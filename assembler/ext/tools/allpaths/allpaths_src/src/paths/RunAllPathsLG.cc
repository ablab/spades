///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "RunAllPathsLG is the assembly pipeline controller for ALLPATHS-LG, the large "
  "genome assembler from the Broad Institute. Although there are many options, "
  "most are experimental and should not be used unless otherwise instructed. "
  "Please use with the default options to avoid running into difficulties.";

/**
   ALLPATHS-LG, the 'Large Genome' version of ALLPATHS

   The various directories in play here:

      ref dir - contains the reference genome we are using, or in the case of a
         unknown genome it just contains the data directories.

      data dir - contains a copy of the genome (if available), and the raw reads
         (unpaired and paired), plus things derived from these data.  These are
         the things that should never change regardless of the RunAllPaths
         options you use.
	 
      run dir - contains everything else, expect the assembly proper.  If you
         need to modify the reads, copy them here first and work on the copy;
         the reads in the data dir must remain the original raw reads. All the
         pre-processing work is done.

      sub dir - contains the actual assembly, built from the pre-processed data
         in the run dir.
*/

#include "MainTools.h"
#include "Basevector.h"
#include "kmers/KmerShape.h"
#include "paths/RunAllPathsCommon.h"
#include "paths/RunAllPathsCommonRules.h"
#include "system/MiscUtil.h"


int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;  
  CommandDoc(DOC);

  // Required Options 

  CommandArgument_String(PRE);
  CommandArgument_String_Doc(DATA_SUBDIR,
    "The data directory name to use.");
  CommandArgument_String(RUN);
  CommandArgument_String_Doc(REFERENCE_NAME, 
    "Organism or reference genome name - used for top level directory name.");

  // Directory Options (DATA and SOURCE) 

  CommandArgument_String_OrDefault_Doc(SUBDIR, "test",
    "directory under REF/DATA/RUN/ASSEMBLIES where the actual assembly will be created");

  // Core Kmer Size
  CommandArgument_Int_OrDefault_Doc(K, 96,
    "core kmer size - only K=96 is actively supported" );

  // Reference Assembly Options

  CommandArgument_String_OrDefault_Valid_Doc(EVALUATION, "BASIC",
    "{NONE,BASIC,STANDARD,FULL,CHEAT}",
    "BASIC - basic evaluations that don't require a reference.\n" 
    "STANDARD - run evaluation modules using a supplied reference.\n"
    "FULL - turn on in-place evaluation in certain assembly modules. "
    "Does not perturb the assembly results, but if used on an existing assembly"
    "it will cause the assembly process to restart part way through.\n"
    "CHEAT - uses the reference to guide the assembly slightly, allowing more detailed"
    "analysis but may cause small but (hopefully) neutral changes to the assembly.");

  // Parallelization Options
  CommandArgument_Int_OrDefault_Doc(MAXPAR, 1,
    "maximum number of concurrent modules to run in pipeline");
  CommandArgument_String_OrDefault_Doc(THREADS, "max", 
    "Number of threads to use in multithreaded modules, unless "
    "overridden by module specific values below. Choices are: "
    "{[1-N], max}");
  CommandArgument_Int_OrDefault_Doc(LR_THREADS, -1,
    "number of parallel threads in LocalizeReadsLG");
  CommandArgument_Int_OrDefault_Doc(CP_THREADS, -1,
    "number of parallel threads in CommonPather");
  CommandArgument_Int_OrDefault_Doc(RDR_THREADS, -1,
    "number of parallel threads in RemoveDodgyReads");
  CommandArgument_Int_OrDefault_Doc(PC_THREADS, -1,
    "number of parallel threads in PreCorrect");
  CommandArgument_Int_OrDefault_Doc(FE_THREADS, -1,
    "number of parallel threads in FindErrors");
  CommandArgument_Int_OrDefault_Doc(FF_THREADS, -1,
    "number of parallel threads in FillFragments");
  CommandArgument_Int_OrDefault_Doc(ECJ_THREADS, -1,
    "number of parallel threads in ErrorCorrectJump");
  CommandArgument_Int_OrDefault_Doc(KP_THREADS, -1,
    "number of parallel threads in FastbToKmerParcels");
  CommandArgument_Int_OrDefault_Doc(CUG_THREADS, -1,
    "number of parallel threads in CloseUnipathGaps");
  CommandArgument_Int_OrDefault_Doc(MN_THREADS, -1,
    "number of parallel threads in MergeNeighborhoods2");
  CommandArgument_Int_OrDefault_Doc(PR_THREADS, -1,
    "number of parallel threads in PathReads");

  // Memory Usage Options
  CommandArgument_LongLong_OrDefault_Doc(MAX_MEMORY_GB,0,
    "try not to use more than this amount of memory - default is to use all memory");

  // Pipeline Options (see also Pipeline Target & Route Options below)

  CommandArgument_Bool_OrDefault_Doc(OVERWRITE, False, 
    "overwrite existing files");
  CommandArgument_Bool_OrDefault_Doc(DRY_RUN, False,
    "build makefile but don't execute it - tests dependances and lets you see "
    "which commands will be excuted without actually running the pipeline");
  CommandArgument_Bool_OrDefault_Doc(VIEW_PIPELINE_AND_QUIT, False,
    "only write out the pipeline description info and quit");
  CommandArgument_Bool_OrDefault_Doc(MEMORY_DIAGNOSTICS, True,
    "Run MemMonitor diagnostics on pipeline modules" );
  CommandArgument_Bool_OrDefault(IGNORE_INITIAL_PARTIAL_OUTPUTS, False);
  CommandArgument_Bool_OrDefault_Doc(IGNORE_METADATA, False,
    "ignores all metadata files - useful if many things are run outside the pipeline");
  CommandArgument_Bool_OrDefault_Doc(USE_LINK_FOR_COPY, True,
    "Use symbolic links instead of copying files" );
  CommandArgument_String_OrDefault_Doc(MAKE_CMD, "make",
    "Make command to use to execute ALLPATHS pipeline makefile");
  CommandArgument_String_OrDefault_Doc(TEST_RULE, "", 
    "Check targets and dependencies are correct for pipeline rules."
    "To test all rules use ALL. To test a single rule use the rule ID or Command name");
  CommandArgument_String_OrDefault_Valid_Doc(TEST_MODE, "standard",
    "{list,prepare,module,standard,full}",
    "list - list targets, dependencies and command\n"
    "prepare - prepare testing directory and list command\n"
    "module - test module, dependencies\n"
    "standard - test module, dependencies, targets\n"
    "full - test module, dependencies, targets, unused dependencies\n");
  CommandArgument_Bool_OrDefault_Doc(TEST_ERASE, True,
    "Erase the testing directory at end of each test - unless TEST_MODE=prepare");

  // Pipeline Targets Options

  CommandArgument_StringSet_OrDefault_Doc(TARGETS, "standard",
    "pseudo targets");
  CommandArgument_String_OrDefault_Doc(TARGETS_REF, "",
    "target files to make in ref_dir");
  CommandArgument_String_OrDefault_Doc(TARGETS_DATA, "",
    "target files to make in data dir");
  CommandArgument_String_OrDefault_Doc(TARGETS_RUN, "",
    "target files to make in run_dir");
  CommandArgument_String_OrDefault_Doc(TARGETS_SUBDIR, "",
    "target files to make in subdir");
  CommandArgument_Bool_OrDefault_Doc(FORCE_TARGETS, False,
    "make these targets even if they seem to be up-to-date.");
  CommandArgument_StringSet_OrDefault_Doc(FORCE_TARGETS_OF, "",
    "make targets of these commands even if they're up-to-date. ");
  CommandArgument_String_OrDefault_Doc(DONT_REBUILD_SUBDIR, "",
    "don't rebuild these files in SUBDIR even if they're out-of-date. "
    "just reuse the existing files. ");
  CommandArgument_StringSet_OrDefault_Doc(DONT_UPDATE_TARGETS_OF, "",
    "don't update targets of these commands -- assume they already exist.");
  CommandArgument_StringSet_OrDefault_Doc(RESTORE_TARGETS_OF, "",
    "restore targets of these commands");
  CommandArgument_Bool_OrDefault_Doc(MAKE, True,
    "only regenerate files if the files they depend on have changed, i.e. act "
    "like that 'make' program. if False, everything is remade afresh." );

  // Pipeline Route Options (which modules to run)

  CommandArgument_Bool_OrDefault_Doc(VALIDATE_INPUTS, True,
    "Perform basic santity checks on the original input data");
  CommandArgument_Bool_OrDefault_Doc(REMOVE_DODGY_READS_FRAG, True,
    "Call RemoveDodgyReads on fragment reads before error correction");
  CommandArgument_Bool_OrDefault_Doc(REMOVE_DODGY_READS_JUMP, True,
    "Call RemoveDodgyReads on jumping reads before error correction");
  CommandArgument_Bool_OrDefault_Doc(PRE_CORRECT, True,
    "call PreCorrect on fragment reads prior to running FindErrors");
  CommandArgument_Bool_OrDefault_Doc(ERROR_CORRECT_FRAGS, True,
    "Error correct fragment reads");
  CommandArgument_Bool_OrDefault_Doc(EXTEND_UNIPATHS, True,
    "Extend fill fragment unipaths using smaller K." );
  CommandArgument_Bool_OrDefault_Doc(CLOSE_UNIPATH_GAPS, True,
    "Run CloseUnipathGaps on unipaths.");
  CommandArgument_Bool_OrDefault_Doc(SIMPLE_CLOSE_GAPS, False,
    "Run SimpleGapCloser on unipaths.");
  CommandArgument_Bool_OrDefault_Doc(LITTLE_HELPS_BIG, True,
    "Run LittleHelpsBig on unipaths.");
  CommandArgument_Bool_OrDefault_Doc(PATCH_UNIPATHS, True,
    "Run UnipathPatcher on unipaths.");
  CommandArgument_Bool_OrDefault_Doc(PATCH_SCAFFOLDS, True,
    "Run PostPatcher on scaffolds.");
  CommandArgument_Bool_OrDefault_Doc(MERGE_UNIPATHS_AS_READS, True,
    "Use unipaths as filled fragment reads in localization." );
  CommandArgument_Bool_OrDefault_Doc(ERROR_CORRECT_LONG_JUMPS, True,
    "Error correct long (20K+) jumping reads?");
  CommandArgument_Bool_OrDefault_Doc(CLOSE_WITH_MIXMERS, False,
    "Close unipath gaps using Mixmers");
  CommandArgument_Bool_OrDefault_Doc(PICK_THE_RIGHT_BRANCH, False,
    "Run PickTheRightBranch");
  CommandArgument_Bool_OrDefault_Doc(UNIPATH_BABY_STATS, False,
    "Run UnipathBabyStats");
  CommandArgument_Bool_OrDefault_Doc(ADD_UNIBASE_JUNCTIONS, False,
    "Add unibase junctions");
  CommandArgument_Bool_OrDefault_Doc(FILTER_PRUNED_READS, True,
    "Filter pruned reads");
  CommandArgument_Bool_OrDefault_Doc(USE_LONG_JUMPS, True,
    "if available, use long (20K+) jump reads for scaffolding. Otherwise ignore.");
  CommandArgument_Bool_OrDefault_Doc(COMPUTE_LONG_JUMP_SEPS, False,
    "Compute long jumping read separations and standard deviations");
  CommandArgument_Bool_OrDefault_Doc(FIX_ASSEMBLY_BASE_ERRORS, False,
    "Fix minor base errors in final assembly");
  CommandArgument_Bool_OrDefault_Doc(FIX_SOME_INDELS, True,
    "Fix minor indels in final assembly");
  CommandArgument_Bool_OrDefault_Doc(KPATCH, True,
    "Patch some gaps in the final assembly by walking through unipath graph");
  CommandArgument_Bool_OrDefault_Doc(REMODEL, True,
    "Recompute gaps after CleanAssembly.");
  CommandArgument_Bool_OrDefault_Doc(FIX_LOCAL, True,
    "Locally reassemble and and attempt to correct errors");
  CommandArgument_Bool_OrDefault_Doc(CLEAN_ASSEMBLY, True,
    "Clean up assembly by removing small contigs and scaffolds");
  CommandArgument_Bool_OrDefault_Doc(CONNECT_SCAFFOLDS, False,
    "Attempt to join scaffolds whose ends overlaps");
  CommandArgument_Bool_OrDefault_Doc(LONG_READ_UNIPATH_PATCH, False,
    "If available, use optional long reads to patch unipaths.");
  CommandArgument_Bool_OrDefault_Doc(LONG_READ_POST_PATCH, True,
    "If available, use optional long reads to patch scaffolded contigs.");
  CommandArgument_Bool_OrDefault_Doc(BIG_MAP, True,
    "Use the BigMap alternative route");

  // ValidateAllPathsInputs Options
  CommandArgument_Bool_OrDefault_Doc(VAPI_WARN_ONLY, False,
    "Validate input data, but continue even if problems are found");

  // FindErrors Options

  CommandArgument_Int_OrDefault_Doc(EC_K, 24,
    "K size to use for error correction in FindErrors");
  CommandArgument_Int_OrDefault_Doc(FE_AUTO_CONFIRM, -1,
    "set AUTO_CONFIRM for FindErrors");
  CommandArgument_Int_OrDefault_Doc(FE_MIN_READSTACK_DEPTH, 5,
    "set MIN_READSTACK_DEPTH for FindErrors");
  CommandArgument_Int_OrDefault_Doc(FE_MIN_BASESTACK_DEPTH, 1,
    "set MIN_BASESTACK_DEPTH for FindErrors");
  CommandArgument_Int_OrDefault_Doc(FE_QUAL_CEIL_RADIUS, 2,
    "set QUAL_CEIL_RADIUS for FindErrors");
  CommandArgument_Bool_OrDefault_Doc(FE_QCR_ALWAYS, True,
    "set QCR_ALWAYS for FindErrors");
  CommandArgument_Int_OrDefault_Doc(FE_MIN_QSS_TO_WIN, 40,
    "set MIN_QSS_TO_WIN for FindErrors");
  CommandArgument_Int_OrDefault_Doc(FE_MIN_QSS_TO_SUPPORT, 60,
    "set MIN_QSS_TO_SUPPORT for FindErrors");
  CommandArgument_Int_OrDefault_Doc(FE_MIN_QSS_TO_CONFIRM, 90,
    "set MIN_QSS_TO_CONFIRM for FindErrors");
  CommandArgument_Int_OrDefault_Doc(FE_MAX_KMER_FREQ_TO_MARK, 2,
    "set MAX_KMER_FREQ_TO_MARK for FindErrors");
  CommandArgument_Int_OrDefault_Doc(FE_NUM_CYCLES, 2,
    "set NUM_CYCLES for FindErrors");
  CommandArgument_Double_OrDefault_Doc(FE_MAX_QSS_RATIO_TO_CORRECT, 0.25,
    "set MAX_QSS_RATIO_TO_CORRECT for FindErrors");
  CommandArgument_Double_OrDefault_Doc(FE_MAX_QSS_RATIO_TO_CORRECT2, 0.25,
    "set MAX_QSS_RATIO_TO_CORRECT2 for FindErrors");
  CommandArgument_Double_OrDefault_Doc(FE_MAX_QSS_RATIO_TO_CONFIRM, 1.0,
    "set MAX_QSS_RATIO_TO_CONFIRM for FindErrors");
  CommandArgument_Bool_OrDefault_Doc(FE_SKIP_UNDER, True,
    "set SKIP_UNDER for FindErrors");
  CommandArgument_Bool_OrDefault_Doc(FE_DO_BRANCHES, False,
    "set DO_BRANCHES for FindErrors");
  CommandArgument_Bool_OrDefault_Doc(FE_USE_KMER_SPECTRUM, False,
    "use kmer spectrum metrics for setting FE parameters");
  CommandArgument_Bool_OrDefault_Doc(HAPLOIDIFY, False,
    "attempt to remove ambiguities/polymorphisms from corrected reads");

  // FillFragment Options

  CommandArgument_Int_OrDefault_Doc(FF_K, 28,
    "kmer size used in fragment filling" );
  CommandArgument_Bool_OrDefault_Doc(FF_UNIQUE_ONLY, False,
    "only accept unique closures (multiple closures not allowed)");
  CommandArgument_Int_OrDefault(FF_MAX_STRETCH, 3.0);
  CommandArgument_Int_OrDefault(FF_MAX_CLIQ, 64000);
  CommandArgument_Int_OrDefault_Doc(FF_MIN_OVERLAP, 28,
    "minimum overlaping bases required by bridging read");
  CommandArgument_String_OrDefault_Doc(FF_EXTRA_FILLS, "",
    "pass to FillFragments as EXTRA_FILLS argument");

  // ErrorCorrectJump Options

  CommandArgument_Int_OrDefault_Doc(ECJ_MIN_MATCH, 20,
    "pass to ErrorCorrectJump as MIN_MATCH");
  CommandArgument_Int_OrDefault_Doc(ECJ_MISMATCH_THRESHOLD, 3,
    "pass to ErrorCorrectJump as MISMATCH_THRESHOLD");
  CommandArgument_Int_OrDefault_Doc(ECJ_MISMATCH_NEIGHBORHOOD, 8,
    "pass to ErrorCorrectJump as MISMATCH_NEIGHBORHOOD");
  CommandArgument_Int_OrDefault_Doc(ECJ_MISMATCH_BACKOFF, 3,
    "pass to ErrorCorrectJump as MISMATCH_BACKOFF");
  CommandArgument_Int_OrDefault_Doc(ECJ_SCORE_DELTA, 20,
    "pass to ErrorCorrectJump as SCORE_DELTA");
  CommandArgument_Int_OrDefault_Doc(ECJ_SCORE_MAX, 100,
    "pass to ErrorCorrectJump as SCORE_MAX");

  // UnibaseCopyNumber Options

  CommandArgument_Bool_OrDefault_Doc(UCN_LOWER, False,
    "pass LOWER to UnibaseCopyNumber program");
  CommandArgument_Bool_OrDefault_Doc(UCN_GAP_REMODEL, True,
    "pass EXP_GAP_REMODEL to UnibaseCopyNumber program");

  // CorrectLongRead Options

  CommandArgument_Int_OrDefault_Doc(CLR_KOUT, 640, "passed to CorrectLongReads");
  CommandArgument_Int_OrDefault_Doc(CLR1_MIN_PATCH1, 200, 
       "passed to first instance of CorrectLongReads");
  CommandArgument_Int_OrDefault_Doc(CLR1_MIN_PATCH2, 200, 
       "passed to first instance of CorrectLongReads");

  // BuildUnipathLinkGraph Options
     
  CommandArgument_Bool_OrDefault_Doc(BULG_TRANSITIVE_FILL_IN, True,
    "for BuildUnipathLinkGraphsLG" );
  CommandArgument_Int_OrDefault_Doc(BULG_MIN_UNIPATH_FOR_TRANSITIVE_JOIN, 100,
    "passed to BuildUnipathLinkGraphsLG");
  CommandArgument_Int_OrDefault_Doc(BULG_LINK_GRAPH_NORMAL_MULT, 10,
    "for BuildUnipathLinkGraphsLG" );

  // LocalizeReadsLG Options

  CommandArgument_String_OrDefault_Doc(LR_PASS_THRU, "",
    "Additional options to pass through to LocalizeReads. Place quotes around args, e.g. "
    "LOCALIZE_READS_PASS_THRU=\"ARG1=value1 ARG2=value2\" etc...");
  CommandArgument_Bool_OrDefault_Doc(LR_EVAL, True,
    "Evaluate neighborhoods within LocalizeReads");
  CommandArgument_String_OrDefault_Doc(LR_LOCAL_DUMP, "",
    "Pass as LOCAL_DUMP to LocalizeReadsLG.");
  CommandArgument_Bool_OrDefault_Doc(LR_READ_CLOUD_EVAL, False,
    "Evaluate read clouds within LocalizeReads");
  CommandArgument_Int_OrDefault_Doc(LR_NHOOD_CN_FACTOR, 1, 
    "Controlling the copy number of unipaths entering localized neighborhood <= "
    "LR_NHOOD_CN_FACTOR * ploidy. If set to 0 then seed cn is the limit.");			     
  CommandArgument_Int_OrDefault_Doc(LR_NHOOD_GLUEPERFECT_SIZE, 0, 
    "If >0 then try to merge duplicate edge sequences of this size or greater");	     
  CommandArgument_Bool_OrDefault_Doc(LR_OVERRIDE_TIMEOUTS, True,
    "pass to LocalizeReadsLG");
  CommandArgument_Bool_OrDefault_Doc(LR_BRIDGE_BETWEEN_CLOUD_UNIPATHS, False,
    "Attempt to walk in the global unipath between cloud unipaths.");
  
  // MergeNeighborhoods Options

  CommandArgument_Int_OrDefault_Doc( MN_MIN_ALIGN_LENGTH, 500, 
    "passed to MergeNeighborhoods");
  CommandArgument_Int_OrDefault_Doc( MN_MIN_OVERLAP, 10000, 
    "passed to MergeNeighborhoods");
  CommandArgument_Int_OrDefault_Doc( MN_MIN_PROPER_OVERLAP, 2000, 
    "passed to MergeNeighborhoods");
  CommandArgument_Int_OrDefault_Doc( MN_MIN_PROPER_OVERLAP_FINAL, 5000, 
    "passed to MergeNeighborhoods");
  
  // MakeScaffoldsLG options

  CommandArgument_Bool_OrDefault( MS_REGAP, True );
  CommandArgument_Bool_OrDefault( MS_REGAP_NEG_GAPS, True );
  CommandArgument_String_OrDefault( MS_PASS_THRU, "" );
  CommandArgument_Bool_OrDefault( MS_FIT_GAPS, True);
  CommandArgument_Bool_OrDefault( MS_FIT_GAPS_EXP, False);
  CommandArgument_Bool_OrDefault( MS_USE_UNIBASES, True);

  // FlattenHKP Options

  CommandArgument_Bool_OrDefault_Doc(FLA_NEW_ALGORITHM, False,
    "use new algorithm for FlattenHKP");

  // PostPatcher options

  CommandArgument_Bool_OrDefault( PP_MOC_STRINGENT, True );

  // FixLocal options

  CommandArgument_Bool_OrDefault( FL_RECYCLE, False );

  // LongReadUnipathPatcher options

  CommandArgument_Bool_OrDefault( LRUP_TERMINAL_ONLY, False );
  CommandArgument_Bool_OrDefault( PUP_NEUTER, False );

  // Assisted patching: use reference to guide contig patching prior to scaffolding
  CommandArgument_Bool_OrDefault(ASSISTED_PATCHING, False);

  // CleanAssembly options

  CommandArgument_Int_OrDefault_Doc(CA_MIN_CONTIG_SIZE_IN, 0,
       "pass as MIN_CONTIG_SIZE_IN to CleanAssembly");
  CommandArgument_Int_OrDefault_Doc(CA_MIN_UNIQUE, 1000,
       "pass as MIN_UNIQUE to CleanAssembly");
  CommandArgument_Int_OrDefault_Doc(CA_MAX_SWALLOW, 500,
       "pass as MAX_SWALLOW to CleanAssembly");
  CommandArgument_Int_OrDefault_Doc(CA_MIN_INDIRECT_GAP_FOR_HYBRID, 800,
       "for hybrid assemblies, pass as MIN_INDIRECT_GAP to CleanAssembly");

  // SubmissionPrep options

  CommandArgument_Bool_OrDefault_Doc(SP_REORDER, False,
    "reorder scaffolds by size after preparing for submission");

  // AssemblyAccuracy options

  CommandArgument_String_OrDefault(AA_REFNAME, "genome" );
  CommandArgument_Int_OrDefault(AA_CHUNK_SIZE, 3000 );

  // ScaffoldAccuracy options

  CommandArgument_String_OrDefault(SA_REFNAME, "genome" );

  // General thresholds

  CommandArgument_Int_OrDefault(MIN_CONTIG, 1000); // note, not used by CleanAssembly
  

  EndCommandArguments;

  // Command Arg validity checks

  ForceAssertSupportedK( K );

  if (RUN == "")  {cout << "You must supply a value for RUN" << endl; exit(1);}
  if (DATA_SUBDIR == "")  {cout << "You must supply a value for DATA_SUBDIR" << endl; exit(1);}

  if ( (PATCH_SCAFFOLDS || FIX_SOME_INDELS) && !PATCH_UNIPATHS )
  {    cout << "If PATCH_SCAFFOLDS or FIX_SOME_INDELS is True, "
	    << "then so must PATCH_UNIPATHS." << endl;
       exit(1);   }

  // Define common directories.
     
  const String base_dir  = PRE;
  String ref_dir   = base_dir + "/" + REFERENCE_NAME;
  while (ref_dir.Contains("//"))
    ref_dir.GlobalReplaceBy( "//", "/" );
  const String data_dir  = ref_dir  + "/" + DATA_SUBDIR;
  const String run_dir   = data_dir + "/" + RUN;
  const String sub_dir   = run_dir  + "/ASSEMBLIES/" + SUBDIR;
  const String log_dir  = ref_dir + "/make_log/" + DATA_SUBDIR + "/" + RUN 
    + "/" + SUBDIR + "/" + (DRY_RUN ? "dry_run" : Date(True /*ISO8601 Format*/) );
  const String run_log_file  = log_dir + "/run.log";
  const String make_dir  = log_dir;

  const String REF = REFERENCE_NAME;
  const String DATA = REF + "/" + DATA_SUBDIR;

  // Define some commonly used strings

  const String KS = ToString( K );
  const String FF_KS = ToString( FF_K );   // Filled Fragment K
  const String dotK = ".k" + KS;
  const String pd = " PRE=" + PRE + " DATA=" + DATA + " ";
  const String pdr  = pd  + " RUN=" + RUN + " ";
  const String pdk  = pd  + " K=" + KS + " ";
  const String pdrk = pdr + " K=" + KS + " ";
  const String pdrs = pdr + " SUBDIR=" + SUBDIR + " ";
  const String pdrsk = pdrk + " SUBDIR=" + SUBDIR + " ";

  // Should we perform full evaluation within assembly modules? (requires reference)
  Bool evalFull = (EVALUATION =="FULL" || EVALUATION == "CHEAT");
  // Should we perform standard evaluation using stand alone modules? (requires reference)
  Bool evalStandard = (EVALUATION =="STANDARD" || evalFull);
  // Should we perform basic evaluation? (reference not required)
  Bool evalBasic = (EVALUATION == "BASIC" || evalStandard);


  // Make sure that if the run directory already exists then OVERWRITE=True
  if ( IsDirectory(run_dir) && !OVERWRITE && !DRY_RUN)
    FatalErr("Run directory already exists. To continue an assembly use OVERWRITE=True");

  // Create log directory where we put makefiles, logs, graphs, etc...
  Mkpath( log_dir );

  // Redirect standard out to log (unless overidden with TEE)
  command.SetOutputRedirection(run_log_file);

  // Save the RunAllPathsLG command for posterity
  Ofstream(commandFile, log_dir + "/command" );
  PrintCommandPretty(commandFile);
  commandFile.close();

  // Warn if using non standard options
  
  if (K != 96)
    cout << "WARNING: Unsupported value of K. Only K=96 is currently supported" << endl << endl;

  if (Member(TARGETS, String("all")))
    cout << "WARNING: Unsupported value of TARGETS. " 
	 << "TARGETS=all is intended for developer use only. " << endl << endl;

  // Display basic machine stats for debugging

  DisplaySystemStatus();

  // Parse threading
  int global_threads =
    ParseThreadArgs(THREADS, LR_THREADS, CP_THREADS, FF_THREADS, ECJ_THREADS,
		    KP_THREADS, CUG_THREADS, RDR_THREADS, FE_THREADS, MN_THREADS, PR_THREADS, 
		    PC_THREADS);
  
  // Create data, run, subdir and tmp directories
  cout << "run dir is \n" << run_dir << "\n" << endl;

  if (!DRY_RUN) {
    Mkpath( sub_dir  );
    ForceAssert( IsDirectory( sub_dir ) );
    // Temporary directories:
    Mkdir777( run_dir + "/tmp" );
    Mkdir777( sub_dir + "/tmp" );
  }

  // Set make

  SetMakeOnlyIfNeeded( MAKE );

  // Check for optional long jumps (used for scaffolding)
  if (USE_LONG_JUMPS && !IsRegularFile( data_dir + "/long_jump_reads_orig.fastb" )) {
    cout << "Unable to find optional long jumping reads (>20kb) for scaffolding." << endl;
    USE_LONG_JUMPS = false;
  }

  // Check for optional long reads (used for patching)
  if ((LONG_READ_POST_PATCH || LONG_READ_UNIPATH_PATCH || BIG_MAP) && 
      !IsRegularFile( data_dir + "/long_reads_orig.fastb" )) {
    cout << "Unable to find optional long reads for patching." << endl;
    LONG_READ_UNIPATH_PATCH = LONG_READ_POST_PATCH = BIG_MAP = false;
  } else {
    cout << "Found optional long reads for patching." << endl;
  }
  if ( BIG_MAP ) FIX_SOME_INDELS = False;
  if ( BIG_MAP ) FIX_LOCAL = False;
  if ( BIG_MAP ) UCN_LOWER = True;

  // Do we have a reference file available?
  Bool haveReference = IsRegularFile(ref_dir + "/genome.fastb");

  // Make sure reference is available if required
  if (evalStandard && !haveReference) {
    FatalErr("Reference genome '" << ref_dir << "/genome.fastb' not found or specified. "
	     "Use EVALUATION=BASIC or NONE if you don't have a reference.");
  }

  // Make sure ploidy file (or values) are available

  if ( !IsRegularFile( data_dir + "/ploidy" ) )
    FatalErr( "Unable to find the ploidy file: " + data_dir + "/ploidy");

  // We need to peek ahead to figure out the ploidy (bit of cheat this)
  const int ploidy = StringOfFile( data_dir + "/ploidy", 1 ).Int( );
  
  if (ploidy != 1 && ploidy != 2)
    FatalErr("Unsupported ploidy value file - must be 1 or 2");  

  // Create or update reference genome.size file
  if (haveReference)
    UpdateRefGenomeSize(ref_dir);

  /////////////////////////////////////////////////////////////////////////////
  //
  //  Setup Make Manager
  //

  // The make object we use to build up the pipeline
  MakeMgr make( make_dir + "/Makefile" );

  make.SetMakeCommand( MAKE_CMD );

  make.SetLnForCp( USE_LINK_FOR_COPY );
  
  make.SetIgnoreMetaData( IGNORE_METADATA );

  make.SetDiagnostics(MEMORY_DIAGNOSTICS);

  make.SetLogDirectory(log_dir);

  // Add directory abbreviations to simplify the makefile
  make.AddDirAbbrev( ref_dir, "REF" );
  make.AddDirAbbrev( data_dir, "DATA" );
  make.AddDirAbbrev( run_dir, "RUN" );
  make.AddDirAbbrev( sub_dir, "SUB" );

  // Add kN abbreviations
  make.AddSubst( ".kN", ".k" + KS );
  make.AddSubst( ".ffkN", ".k" + FF_KS );

  // Add pre-existing files
  make.AddPreExistingFiles(FilesIn( data_dir, "frag_reads_orig.{fastb,qualb,pairs}",
				    "jump_reads_orig.{fastb,qualb,pairs}","ploidy" ) );

  if (haveReference)
    make.AddPreExistingFiles(FilesIn( ref_dir, "genome.{fasta,fastb,size}" ) );

  make.AddPreExistingFiles( FilesIn( sub_dir, DONT_REBUILD_SUBDIR ) );

  if (USE_LONG_JUMPS)
    make.AddPreExistingFiles(FilesIn( data_dir, "long_jump_reads_orig.{fastb,qualb,pairs}") );

  if (LONG_READ_UNIPATH_PATCH || LONG_READ_POST_PATCH || BIG_MAP)
    make.AddPreExistingFiles(FilesIn( data_dir, "long_reads_orig.fastb") );

  if (AA_REFNAME != "genome")
    make.AddPreExistingFiles(FilesIn( ref_dir, AA_REFNAME + ".fasta" ) );

  if (UNIPATH_BABY_STATS)
    // Require extended genome
    make.AddPreExistingFiles(FilesIn( ref_dir, "genome_extended.fasta" ) );
     
  make.AddDontUpdateTargetsOf( DONT_UPDATE_TARGETS_OF );
  make.AddRestoreTargetsOf( RESTORE_TARGETS_OF );
     
  make.AddIgnoreChangesToCommandArgs( "{PRE,NUM_THREADS,MAX_MEMORY_GB}" );

  make.SetDryRun( DRY_RUN );
  make.SetIgnoreInitialPartialOutputs( IGNORE_INITIAL_PARTIAL_OUTPUTS );

  /////////////////////////////////////////////////////////////////////////////
  //
  //  Make Manager Rule Definitions
  //

  make.PrepareForRules();

  String prepareGenomeTargets = "";

  if ( haveReference ) {
    AddRule_MakeLookupTable(make, "genome.fastb", ref_dir,		     
			    "Making lookup table for reference genome" );

    prepareGenomeTargets =
      "genome.{paths{,_rc,db_big}.kN,unipaths{,db_big,.true_count,.placements}.kN,"
      "unipath_adjgraph.kN,unibases.kN{,.fastb,.lookup}}";

// MakeDepend: dependency PrepareGenome
    make.AddRule( FilesIn( ref_dir, prepareGenomeTargets),
		  FilesIn( ref_dir, "genome.fastb" ),
		  "PrepareGenome PRE=" + PRE + 
		  " DATA=" + REF +
		  " K=" + KS +
		  " CANONICALIZE=True" +
		  " NUM_THREADS=" + ToString(CP_THREADS) ,
		     
		  "Computing many useful things about the reference (e.g. true "
		  "unipaths and their true copy numbers), "
		  "for evaluating how well we approximate these things from the "
		  "reads.  Note that the kmer numbering will differ from any kmer "
		  "numbering in the reads, "
		  "since this was computed without regard to any read set!" );
    
    // Copy genome files into the data directory
    make.AddCpRule( data_dir, ref_dir, "genome.{fastb,fasta,lookup}" );
    make.AddCpRule( data_dir, ref_dir, prepareGenomeTargets );

    if (AA_REFNAME != "genome")
      make.AddCpRule( data_dir, ref_dir, AA_REFNAME +".fasta" );

  }

  // Validate input data

// MakeDepend: dependency ValidateAllPathsInputs
  make.AddOutputRule( FilesIn( data_dir, "ValidateAllPathsInputs.report",
			       "ValidateAllPathsInputs_core_stats.out"),
		      FilesIn( data_dir, "frag_reads_orig.{fastb,pairs,qualb}",
			       "jump_reads_orig.{fastb,pairs,qualb}",
			       (USE_LONG_JUMPS ? "long_jump_reads_orig.{fastb,pairs,qualb}" : "")),
		      "ValidateAllPathsInputs " + pdk +
		      " LONG_JUMPS=" + ToStringBool(USE_LONG_JUMPS) +
		      " WARN_ONLY=" + ToStringBool(VAPI_WARN_ONLY) ,
		      
		      "Performs basic validation on inputs to ALLPATHS." );

  // copy reads (and aux read files) to run directory
  make.AddCpRule( run_dir, data_dir, "frag_reads_orig.{fastb,qualb,pairs}" );
  make.AddCpRule( run_dir, data_dir, "jump_reads_orig.{fastb,qualb,pairs}" );
  make.AddCpRule( run_dir, data_dir, "ploidy" );

  if (LONG_READ_UNIPATH_PATCH | LONG_READ_POST_PATCH || BIG_MAP)
    make.AddCpRule( run_dir, data_dir, "long_reads_orig.fastb" );




  ///////////////////////////////////////////////////////////////////////////
  // Add long jump if they are available
  //


  String long_jump_reads_current = "long_jump_reads_orig";

  if (USE_LONG_JUMPS) {
    make.AddCpRule( run_dir, data_dir, "long_jump_reads_orig.{fastb,qualb,pairs}" );

// MakeDepend: dependency RemoveDodgyReads
    if ( REMOVE_DODGY_READS_JUMP ) {
      make.AddRule( FilesIn( run_dir, "long_jump_reads_filt.{fastb,qualb,pairs,readtrack}" ),
		    FilesIn( run_dir, "long_jump_reads_orig.{fastb,qualb,pairs}" ),
		  "RemoveDodgyReads"
                    " IN_HEAD=" + run_dir + "/long_jump_reads_orig"
                    " OUT_HEAD=" + run_dir + "/long_jump_reads_filt"
                    " REMOVE_DUPLICATES=True"
		    " RC=False"  
		    " NUM_THREADS=" + ToString(RDR_THREADS),
		    "Remove reads that are poly-A or duplicate molecules." );
    }  else {
      make.AddCpRule( run_dir, "long_jump_reads_filt",
		      run_dir, "long_jump_reads_orig",
		      ".{fastb,qualb,pairs}");
    }

    long_jump_reads_current = "long_jump_reads_filt";
  }

// MakeDepend: dependency RemoveDodgyReads
  if ( REMOVE_DODGY_READS_FRAG ) {
    make.AddRule( FilesIn( run_dir, "frag_reads_filt.{fastb,qualb,pairs,readtrack}" ),
		  FilesIn( run_dir, "frag_reads_orig.{fastb,qualb,pairs}" ),
		  "RemoveDodgyReads"
                  " IN_HEAD=" + run_dir + "/frag_reads_orig"
                  " OUT_HEAD=" + run_dir + "/frag_reads_filt"
		  " REMOVE_DUPLICATES=False"
                  " RC=False"  // the reads have NOT been RC'ed
		  " NUM_THREADS=" + ToString(RDR_THREADS),
		  "Remove reads that are poly-A." );
  }
  else {
    make.AddCpRule( run_dir, "frag_reads_filt",
		    run_dir, "frag_reads_orig",
		    ".{fastb,qualb,pairs}");
  }
    
  
  // This block is for any up-front read filtering, QC, etc.  Takes
  // *_reads_orig and produces *_reads_filt. --bruce

  if (REMOVE_DODGY_READS_JUMP) {
    
    make.AddRule( FilesIn( run_dir, "jump_reads_filt.{fastb,qualb,pairs,readtrack}" ),
		  FilesIn( run_dir, "jump_reads_orig.{fastb,qualb,pairs}" ),
		  "RemoveDodgyReads"
                  " IN_HEAD=" + run_dir + "/jump_reads_orig" 
                  " OUT_HEAD=" + run_dir + "/jump_reads_filt"
		  " REMOVE_DUPLICATES=True"
                  " RC=True"  // the reads have been RC'ed
		  " NUM_THREADS=" + ToString(RDR_THREADS),
		  "Remove reads that are poly-A or duplicate molecules." );

  } else {
    make.AddCpRule( run_dir, "jump_reads_filt",
		    run_dir, "jump_reads_orig",
		    ".{fastb,qualb,pairs}");
  }


  if (ERROR_CORRECT_FRAGS) {
    
    const unsigned K_ODD = EC_K & 1 ? EC_K : EC_K + 1;

// MakeDepend: dependency FindErrors
    make.AddRule( FilesIn( run_dir, 
                           "frag_reads_filt." + ToString(K_ODD) + "mer.kspec",
                           "frag_reads_edit.{fastb,qualb}" ),
		  FilesIn( run_dir, 
                           "ploidy",
                           "frag_reads_filt.{fastb,qualb}" ),
                  String("FindErrors")
                  + " K=" + ToString(EC_K)
		  + " HEAD_IN="  + run_dir + "/frag_reads_filt"
		  + " HEAD_OUT=" + run_dir + "/frag_reads_edit"
                  + " PLOIDY_FILE=" + run_dir + "/ploidy"
                  + ARG(PRE_CORRECT, PRE_CORRECT)
                  + ARG(AUTO_CONFIRM, FE_AUTO_CONFIRM)
                  + ARG(MIN_BASESTACK_DEPTH, FE_MIN_BASESTACK_DEPTH)
                  + ARG(MIN_READSTACK_DEPTH, FE_MIN_READSTACK_DEPTH)
                  + ARG(QUAL_CEIL_RADIUS, FE_QUAL_CEIL_RADIUS)
                  + ARG(QCR_ALWAYS, FE_QCR_ALWAYS)
                  + ARG(MIN_QSS_TO_WIN, FE_MIN_QSS_TO_WIN)
                  + ARG(MIN_QSS_TO_SUPPORT, FE_MIN_QSS_TO_SUPPORT)
                  + ARG(MIN_QSS_TO_CONFIRM, FE_MIN_QSS_TO_CONFIRM)
                  + ARG(MAX_QSS_RATIO_TO_CORRECT,  FE_MAX_QSS_RATIO_TO_CORRECT)
                  + ARG(MAX_QSS_RATIO_TO_CORRECT2, FE_MAX_QSS_RATIO_TO_CORRECT2)
                  + ARG(MAX_QSS_RATIO_TO_CONFIRM, FE_MAX_QSS_RATIO_TO_CONFIRM)
                  + ARG(NUM_CYCLES, FE_NUM_CYCLES)
                  + ARG(SKIP_UNDER, FE_SKIP_UNDER)
                  + ARG(DO_BRANCHES, FE_DO_BRANCHES)
		  + ARG(USE_KMER_SPECTRUM, FE_USE_KMER_SPECTRUM)
		  + " NUM_THREADS=" + ToString(FE_THREADS)
                  + " MAX_MEMORY_GB=" + ToString(MAX_MEMORY_GB),

		  "Error correct reads");

    make.AddCpRule( run_dir + "/frag_reads_edit.pairs", run_dir + "/frag_reads_filt.pairs" );



// MakeDepend: dependency CleanCorrectedReads
    make.AddRule( FilesIn( run_dir, 
                           "frag_reads_corr.{fastb,qualb,pairs,readtrack," + ToString(K_ODD) + "mer.kspec}",
                           "frag_reads_edit." + ToString(EC_K) + "mer.kspec",
                           (HAPLOIDIFY ? 
                            "frag_reads_corr.{diploid.25mer.kspec,k25.poly_stats,poly.k25.{in,out}.fasta}" : "") ),
		  FilesIn( run_dir, 
                           "ploidy",
                           "frag_reads_edit.{fastb,qualb,pairs}",
                           (FE_USE_KMER_SPECTRUM ? "frag_reads_filt." + ToString(K_ODD) + "mer.kspec" : "") ),
		  String("CleanCorrectedReads ")
                  + " K=" + ToString(EC_K)
                  + " PLOIDY_FILE=" + run_dir + "/ploidy"
		  + " DELETE=True"
		  + " HEAD_IN=" + run_dir + "/frag_reads_edit"
		  + " HEAD_OUT=" + run_dir + "/frag_reads_corr"
                  + ARG(MAX_KMER_FREQ_TO_MARK, FE_MAX_KMER_FREQ_TO_MARK)
                  + ARG(HAPLOIDIFY, HAPLOIDIFY)
		  + (FE_USE_KMER_SPECTRUM ? 
                     " HEAD_KMER_SPECTRUM_IN=" + run_dir + "/frag_reads_filt" : "")
		  + " NUM_THREADS=" + ToString(FE_THREADS)
                  + " MAX_MEMORY_GB=" + ToString(MAX_MEMORY_GB),

		  "Remove reads with low frequency kmers and maybe haploidify");



// MakeDepend: dependency EvaluateCorrectedPairs
    if ( evalStandard )
      make.AddOutputRule( FilesIn( run_dir, "frag_reads_corr.EvaluateCorrectedPairs.out" ),
			  JoinVecs( FilesIn( run_dir, "frag_reads_corr.{fastb,pairs}",
					     "frag_reads_orig.{fastb,pairs}"),
				    FilesIn( data_dir, "genome.{fastb,lookup}") ),
			  "EvaluateCorrectedPairs " + pdr
			  + " READS_IN=frag_reads_corr READS_ORIG=frag_reads_filt",
			  
			  "Evaluate error corrected paired reads using reference." );

  } else {
    // No Error Correction - use uncorrected reads

    make.AddCpRule( run_dir, "frag_reads_corr", 
		    run_dir, "frag_reads_filt",
		    ".{fastb,qualb,pairs}");
    make.AddCpRule( run_dir, "frag_reads_edit", 
		    run_dir, "frag_reads_filt",
		    ".{fastb}");
  }

  String correctedFragReads = "frag_reads_corr";
  String correctedFragReadsDB = correctedFragReads + ".k" + ToString(FF_KS) +
    ".{ug.{info,seq,edges,dict},pc.{info,seq,paths,ev}}";
// MakeDepend: dependency PathReads
  make.AddRule( FilesIn(run_dir,correctedFragReadsDB),
                FilesIn(run_dir,correctedFragReads+".fastb"),

                "PathReads " + pdr +
                " READS_IN=" + correctedFragReads +
                " NUM_THREADS=" + ToString(PR_THREADS) +
                " MAX_MEMORY_GB=" + ToString(MAX_MEMORY_GB),
                "Path corrected fragment reads");

  ///////////////////////////////////////////////////////////////////////////
  // Create filled reads from fragment libraries
  // Compute frag read pair separations

// MakeDepend: dependency FillFragments
  make.AddRule( FilesIn( run_dir, "filled_reads.{fastb,pairs,info}",
			 correctedFragReads + "_cpd.pairs"),
		FilesIn( run_dir, correctedFragReadsDB, 
			 "frag_reads_corr.pairs"),
		
		"FillFragments " + pdr +
		" READS_IN=" + correctedFragReads +
		" READS_OUT=filled_reads"
		" PAIRS_OUT=" + correctedFragReads + "_cpd"
		" NUM_THREADS=" + ToString(FF_THREADS) +
		" MAX_STRETCH=" + ToString(FF_MAX_STRETCH) +
		" MAX_CLIQ=" + ToString(FF_MAX_CLIQ) +
		" MIN_OVERLAP=" + ToString(FF_MIN_OVERLAP) +
                ARG(EXTRA_FILLS, FF_EXTRA_FILLS) +
                " UNIQUE_ONLY=" + ToStringBool(FF_UNIQUE_ONLY),

		
		"Use corrected fragment reads to create filled fragment reads" );

  String filledReads = "filled_reads";

  make.AddCpRule( run_dir, correctedFragReads + "_cpd",
		  run_dir, correctedFragReads,
		  ".{fastb,qualb}");

  correctedFragReads += "_cpd";

  // copy new pairing info back to uncorrected but filtered frag reads

// MakeDepend: dependency ReplacePairsStats
  make.AddRule( FilesIn( run_dir, "frag_reads_filt_cpd.pairs" ),
	        FilesIn( run_dir, "frag_reads_filt.pairs", correctedFragReads + ".pairs" ),
		"ReplacePairsStats " 
		" PAIRS_IN=" + run_dir + "/frag_reads_filt"
		" PAIRS_OUT=" + run_dir + "/frag_reads_filt_cpd"
		" STATS_IN=" + run_dir + "/" + correctedFragReads,
	
		"Update filtered read stats with new values" );

  make.AddCpRule( run_dir, "frag_reads_filt_cpd",
		  run_dir, "frag_reads_filt",
		  ".{fastb,qualb}");


  ///////////////////////////////////////////////////////////////////////////
  // Path filled fragments
  //

  AddRule_CommonPather(make, CP_THREADS, K, run_dir,
		       run_dir, filledReads + ".fastb");

// MakeDepend: dependency MakeRcDb
  make.AddRule( FilesIn( run_dir, filledReads + ".paths{_rc,db}.kN" ),
		FilesIn( run_dir, filledReads + ".paths.kN" ),
		"MakeRcDb " + pdrk + 
		" READS=" + filledReads );
  

  ///////////////////////////////////////////////////////////////////////////
  // Build filled fragments unipaths
  //

// MakeDepend: dependency Unipather
  make.AddRule( FilesIn( run_dir, filledReads + ".{unipaths{,db},unibases,unipath_adjgraph}.kN" ),
		FilesIn( run_dir, filledReads + ".{paths{,_rc,db}.kN,fastb}"),
		"Unipather " + pdrk +
		" READS=" + filledReads + 
		" UNIPATHS=unipaths UNIBASES=unibases",
		
		"Creating unipaths from filled reads" );


  // Final Unibases and filled reads filenames
  String unibasesName = filledReads;

  String current_unipaths = filledReads;

  ///////////////////////////////////////////////////////////////////////////
  // Extend filled fragment unipaths using lower K and fragment reads pairs
  // and optionally also semi-corrected jumping read pairs
  //

  String final_uni_for_UnipathPatcher;
  if (EXTEND_UNIPATHS) {

    // Use error corrected and separation computed fragment reads from now on
    String fragReads = "frag_reads_corr_cpd";

    // Precorrect unipaths using long reads.

  if ( BIG_MAP )
  {
// MakeDepend: dependency FindUnipathGaps
	make.AddRule( FilesIn( run_dir, current_unipaths + ".unibases.k" 
			       + ToString(K) + ".predicted_gaps.txt" ),
		      FilesIn( run_dir, current_unipaths+ ".unibases.kN",
			     "jump_reads_filt.{pairs,fastb}" ),
		      "FindUnipathGaps " + pdr +
		      ARG(K1, K) +
		      " NUM_THREADS=" + ToString(global_threads) +
		      " JUMP_READS_IN=jump_reads_filt" +
                      " MIN_KMERS=100"
		      " HEAD=" + current_unipaths,
		  
		      "Compute gaps between unipaths" );

// MakeDepend: dependency UnibaseCopyNumber3
  make.AddRule( FilesIn( run_dir, 
			 current_unipaths + ".unipaths.cn_raw.kN" ),
		FilesIn( run_dir, current_unipaths + ".unibases.kN", 
                         correctedFragReads + ".fastb", "ploidy"),
		
		"UnibaseCopyNumber3 " + pdrk +
		" NUM_THREADS=" + ToString(global_threads) +
		" READS=" + current_unipaths +
		" READS_EC=\"{" + correctedFragReads + "}\"" 
		" ERR_RATE=0.0 "
		" LOWER=False "
		" WRITE_GAP=True"
		" EXP_GAP_REMODEL=True",
		
		"Estimating the copy number of unibases (and the corresponding "
		"unipaths), based on how many k-mers pile on each unibase." );   

// MakeDepend: dependency CorrectLongReads

      make.AddRule( FilesIn( run_dir, current_unipaths 
                             + ".CLR_modified.{unibases}.kN" ),
		    FilesIn( run_dir, "frag_reads_edit.{fastb,qualb,pairs}",
                             "long_reads_orig.fastb",
			     current_unipaths 
                                  + ".{unibases.kN,unibases.kN.predicted_gaps.txt,unipaths.cn_raw.kN}" ),
		    "CorrectLongReads " + pdrk +
		    " NUM_THREADS=" + ToString(CUG_THREADS) +
		    " IN_HEAD=filled_reads" +
                    " MIN_TO_PATCH=10" +
                    ARG(MIN_PATCH1, CLR1_MIN_PATCH1) +
                    ARG(MIN_PATCH2, CLR1_MIN_PATCH2) +
                    " CLUSTER_ALIGNS_NEW_CLEAN=True PATCHES=True"
                    " WRITE_MODIFIED_UNIBASES=True",

		    "Carry out initial correction of unipaths using long reads");

      current_unipaths += ".CLR_modified";
    }
    
    // Make unipaths bigger by closing gaps

    if ( CLOSE_UNIPATH_GAPS ) {

      if ( BIG_MAP ) {
// MakeDepend: dependency FindUnipathGaps
	make.AddRule( FilesIn( run_dir, current_unipaths + ".unibases.k" 
			       + ToString(K) + ".predicted_gaps.txt" ),
		      FilesIn( run_dir, current_unipaths+ ".unibases.kN",
			     "jump_reads_filt.{pairs,fastb}" ),
		      "FindUnipathGaps " + pdr +
		      ARG(K1, K) +
		      " NUM_THREADS=" + ToString(global_threads) +
		      " JUMP_READS_IN=jump_reads_filt" +
		      " HEAD=" + current_unipaths,
		  
		      "Compute gaps between unipaths" );
      }

// MakeDepend: dependency CloseUnipathGaps
      make.AddRule( FilesIn( run_dir, current_unipaths + ".gap_closed.{unipaths,unipathsdb,unibases,unipath_adjgraph}.kN",
			     current_unipaths + ".gap_closed.{CloseUnipathGaps.log,unibases.plus_junctions.kN}"),
		    FilesIn( run_dir, fragReads + ".{fastb,qualb,pairs}",
			     current_unipaths + ".{unibases.kN"
			     + ( BIG_MAP ?  ",unibases.kN.predicted_gaps.txt" : "" ) 
			     + "}" ),
		    "CloseUnipathGaps "
		    " NUM_THREADS=" + ToString(CUG_THREADS) +
		    " DIR=" + run_dir +
		    " IN_HEAD=" + fragReads +
		    " UNIBASES=" + current_unipaths+ ".unibases.kN"
		    " UNIBASES_K=" + KS +
                    " USE_LINKS=" + ( BIG_MAP ? "True" : "False" ) +
		    " OUT_HEAD=" + current_unipaths + ".gap_closed"
		    " WORKDIR=tmp",
		    
		    "Try extending reads to close gaps between unipaths");

      current_unipaths += ".gap_closed";
    }

    // Make unipaths bigger by closing gaps without using pairing.

// MakeDepend: dependency SimpleGapCloser
    if ( SIMPLE_CLOSE_GAPS ) {
      make.AddRule( FilesIn( run_dir, current_unipaths + ".simple_closed"
			     ".{paths,paths_rc,unipaths,pathsdb,unipathsdb,unibases}.kN" ),
		    FilesIn( run_dir, "frag_reads_filt.{fastb,qualb,pairs}", 
                          "jump_reads_filt.{fastb,qualb,pairs}", 
			     current_unipaths+ ".unibases.kN" ),
		    "SimpleGapCloser " + pdrk +
		    " NUM_THREADS=" + ToString(CUG_THREADS) +
		    " HEAD_IN=" + current_unipaths,
		    " HEAD_OUT=" + current_unipaths + ".simple_closed",
		    
		    "Try directly closing gaps between unipaths");

      current_unipaths += ".simple_closed";
    }

    if (BIG_MAP) {

    AddRule_CommonPather(make, CP_THREADS, 40, run_dir,
			 run_dir, fragReads + ".fastb");
      
// MakeDepend: dependency MakeRcDb
    make.AddRule( FilesIn( run_dir, fragReads + ".paths{_rc,db}.k40" ),
		  FilesIn( run_dir, fragReads + ".paths.k40" ),
		  "MakeRcDb " + pdr + "K=40" +
		  " READS=" + fragReads );
      
// MakeDepend: dependency Unipather
    make.AddRule( FilesIn( run_dir, fragReads + ".{unipaths{,db},unibases}.k40" ),
		  FilesIn( run_dir, fragReads + ".{paths{,_rc,db}.k40,fastb}"),
		  "Unipather " + pdr +
		  " K=40"
		  " BUILD_UNIGRAPH=False"
		  " READS=" + fragReads +
		  " UNIPATHS=unipaths UNIBASES=unibases",
		  
		  "Creating K40 unipaths from edited frag reads" );
    } else {
// MakeDepend: dependency FastbToUnibases
    make.AddRule( FilesIn( run_dir, fragReads + ".unibases.k40" ),
		  FilesIn( run_dir, fragReads + ".fastb"),
		  "FastbToUnibases "
		  " FASTB_IN=" + run_dir + "/" + fragReads + ".fastb" 
		  " NUM_THREADS=" + ToString(global_threads),
		  
		  "Creating K40 unipaths from edited frag reads" );      
    }
      
    if (LITTLE_HELPS_BIG) {
      if (BIG_MAP) {
// MakeDepend: dependency LittleHelpsBig
    make.AddRule( FilesIn( run_dir, "extended40.{unipaths,unipathsdb,unibases,unipath_adjgraph}.kN"),
		  FilesIn( run_dir, fragReads + ".unibases.k40",
			   current_unipaths + ".{unipaths,unibases,unipath_adjgraph}.kN"),
		  "LittleHelpsBig " 
		  " IN_HEAD_SM=" + run_dir + "/" + fragReads +
		  " K_SM=40"
                  " IN_HEAD_LG=" + run_dir + "/" + current_unipaths +
		  " K_LG=" + KS + 
                  " OUT_HEAD=" + run_dir + "/extended40"
		  " NUM_THREADS=" + ToString(CP_THREADS),
		  
		  "Extend unipaths with K=40 unipaths");

      } else {
// MakeDepend: dependency LittleHelpsBig
    make.AddRule( FilesIn( run_dir, "extended40.part1.fastb"),
		  FilesIn( run_dir, fragReads + ".unibases.k40",
			   current_unipaths + ".{unipaths,unibases,unipath_adjgraph}.kN"),
		  "LittleHelpsBig " 
		  " REPATH=False"
		  " IN_HEAD_SM=" + run_dir + "/" + fragReads +
		  " K_SM=40"
                  " IN_HEAD_LG=" + run_dir + "/" + current_unipaths +
		  " K_LG=" + KS + 
                  " OUT_HEAD=" + run_dir + "/extended40"
		  " NUM_THREADS=" + ToString(CP_THREADS),
		  
		  "Extend unipaths with K=40 unipaths");

// MakeDepend: dependency LittleHelpsBigPart2
    make.AddRule( FilesIn( run_dir, "extended40.{unipaths,unipathsdb,unibases,unipath_adjgraph}.kN"),
		  FilesIn( run_dir, current_unipaths + ".unipaths.kN",
			   "extended40.part1.fastb"),
		  "LittleHelpsBigPart2 " 
                  " IN_HEAD_LG=" + run_dir + "/" + current_unipaths +
		  " K_LG=" + KS + 
                  " OUT_HEAD=" + run_dir + "/extended40"
		  " NUM_THREADS=" + ToString(CP_THREADS),
		  
		  "Repathing outside LHB to avoid SIGBUS problem - tmp solution");
      }

    } else {
      make.AddCpRule( run_dir, "extended40", run_dir, current_unipaths,
		      ".{unipaths,unipathsdb,unibases,unipath_adjgraph}.kN" );
    }

    // resetting unipath name
    current_unipaths = "extended40";


// MakeDepend: dependency ShaveUnipathGraph
    make.AddRule( FilesIn( run_dir, current_unipaths +
			   ".shaved.{unipaths,unipathsdb,unibases,unipath_adjgraph}.kN"),
		  FilesIn( run_dir, current_unipaths +
			   ".{unipaths,unipathsdb,unibases,unipath_adjgraph}.kN"),
		  "ShaveUnipathGraph "
		  " DIR=" + run_dir +
		  " IN_HEAD=" + current_unipaths +
		  " OUT_HEAD=" + current_unipaths + ".shaved"
		  " K=" + KS +
		  " NUM_THREADS=" + ToString(CP_THREADS),
		  
		  "Shave extended unipaths (K40)");
    
    current_unipaths += ".shaved";

    ///////////////////////////////////////////////////////////////////////////
    // Compute jumping read separations
    //

// MakeDepend: dependency SamplePairedReadStats
    make.AddRule( FilesIn( run_dir, "jump_reads_filt_cpd.pairs", 
			   "jump_reads_filt.{outies,sample_stats}" ),
		  JoinVecs(  FilesIn( run_dir, "jump_reads_filt.{fastb,qualb,pairs}" ),
			     FilesIn( run_dir, current_unipaths + ".unibases.kN" ) ),  
		  "SamplePairedReadStats " + pdr +
		  " NUM_THREADS=" + ToString(global_threads) +
		  " UNIBASES_K=" + KS +
		  " UNIBASES=" + current_unipaths + ".unibases.k" + KS +
		  " READS_HEAD=jump_reads_filt" 
		  " OUT_HEAD=jump_reads_filt_cpd"
		  " WRITE_SEPS=False",
		      
		  "Use filled fragment unipaths to compute read pair separations" );

    make.AddCpRule( run_dir, "jump_reads_filt_cpd", 
		    run_dir, "jump_reads_filt", 
		    ".{fastb,qualb,outies}" );
    
    if (USE_LONG_JUMPS && COMPUTE_LONG_JUMP_SEPS) { 

      make.AddRule( FilesIn( run_dir, "long_jump_reads_filt_cpd.pairs",
			     "long_jump_reads_filt.{outies,sample_stats}" ),
		    JoinVecs(  FilesIn( run_dir, "long_jump_reads_filt.{fastb,qualb,pairs}" ),
			       FilesIn( run_dir, current_unipaths + ".unibases.kN" ) ),  
		    "SamplePairedReadStats " + pdr +
		    " NUM_THREADS=" + ToString(global_threads) +
		    " UNIBASES_K=" + KS +
		    " UNIBASES=" + current_unipaths + ".unibases.k" + KS +
		    " READS_HEAD=long_jump_reads_filt" 
		    " OUT_HEAD=long_jump_reads_filt_cpd"
		    " WRITE_SEPS=False",
		    
		    "Use filled fragment unipaths to compute read pair separations" );

      make.AddCpRule( run_dir, "long_jump_reads_filt_cpd", 
		      run_dir, "long_jump_reads_filt", 
		      ".{fastb,qualb,outies}" );

      long_jump_reads_current = "long_jump_reads_filt_cpd";
    }

  
    final_uni_for_UnipathPatcher = current_unipaths;
    if (PATCH_UNIPATHS) {
// MakeDepend: dependency UnipathPatcher
      make.AddRule( FilesIn( run_dir, current_unipaths + ".patched." + "unibases" + ".kN",
			     "unipath_patch/UnipathPatcher.{ALIGNS,JALIGNS,SEGS,JSEGS,JOINDATA}",
			     "unipath_patch/extenders_{b,q}_frag_file",
			     "unipath_patch/extenders_{b,q}_jump_file"),
		    FilesIn( run_dir, current_unipaths + ".unibases.kN",
			     "frag_reads_filt_cpd.{fastb,qualb,pairs}",
			     "jump_reads_filt_cpd.{fastb,qualb,pairs,outies}",
			     "frag_reads_edit.fastb", "ploidy" ),
		     "UnipathPatcher " + pdrk
		    + ARG(NUM_THREADS, global_threads)
		    + " IN_HEAD=" + current_unipaths
		    + " OUT_HEAD=" + current_unipaths + ".patched"
                    + ( (CLOSE_WITH_MIXMERS || LONG_READ_UNIPATH_PATCH || BIG_MAP) ? 
			" BONA_FIDE_UNIPATHS=True" : "" )
		    + " MAX_MEMORY_GB=" + ToString(MAX_MEMORY_GB)
		    + ARG(FRAG_READS, "frag_reads_filt_cpd")
		    + ARG(JUMP_READS, "jump_reads_filt_cpd") ,
		    
		    "Try to patched unipaths" );

      current_unipaths += ".patched";
    }

    if (LONG_READ_UNIPATH_PATCH) {
// MakeDepend: dependency LongReadUnipathPatcher
      make.AddRule( FilesIn( run_dir, current_unipaths + ".longread." + "unibases" + ".kN" ),
		    FilesIn( run_dir, current_unipaths + ".unibases.kN",
			     "long_reads_orig.fastb"),
		    "LongReadUnipathPatcher " + pdrk
		    + ARG(NUM_THREADS, global_threads)
		    + " IN_HEAD=" + current_unipaths
		    + " OUT_HEAD=" + current_unipaths + ".longread"
		    + " READS=long_reads_orig"
		    + " PATCH_ONLY=True"
                    + ARG(TERMINAL_ONLY, LRUP_TERMINAL_ONLY)
                    + ARG(NEUTER, PUP_NEUTER),

		    "Try to patched unipaths with longread data" );

      current_unipaths += ".longread";
    }

  String correctedJumpReads = "jump_reads_ec";


    ///////////////////////////////////////////////////////////////////////////
    // Experimental Unipath modules
    //
    // - warning, the rules are incomplete
    

    if (CLOSE_WITH_MIXMERS) {
// MakeDepend: dependency Mixmer
      make.AddRule( FilesIn( run_dir, current_unipaths + ".mixmered" + 
			     ".{unibases,paths,paths_rc,unipaths,pathsdb,unipathsdb}.kN"),
		    FilesIn( run_dir, current_unipaths + ".unibases.kN",
			     "frag_reads_filt.{pairs,fastb,qualb}",
			     "frag_reads_corr.{pairs,fastb}" ),
		    "Mixmer " + pdrk
		    + ARG(NUM_THREADS, global_threads)
		    + ARG(UNIBASES_IN, current_unipaths)
		    + ARG(UNIBASES_OUT, current_unipaths + ".mixmered"),
		    
		    "Try to close gaps with Mixmers" );
      
      current_unipaths += ".mixmered";
    }

    if (PICK_THE_RIGHT_BRANCH) {
// MakeDepend: dependency PickTheRightBranch
      make.AddRule( FilesIn( run_dir, current_unipaths + ".picked" + 
			     ".{unibases,paths,paths_rc,unipaths,pathsdb,unipathsdb}.kN"),
		    FilesIn( run_dir, current_unipaths + ".unibases.kN",
			     "frag_reads_edit.{pairs,fastb,qualb}" ),
		    "PickTheRightBranch " + pdrk
		    + ARG(NUM_THREADS, global_threads)
		    + ARG(UNIBASES_IN, current_unipaths)
		    + ARG(UNIBASES_OUT, current_unipaths + ".picked"),
		    
		    "Clean up unipaths by picking correct branch at junctions" );
      
      current_unipaths += ".picked";
    }

    ///////////////////////////////////////////////////////////////////////////
    // Compute additional unipath metrics
    //
    // - warning, the rules are incomplete.
    // - not fully integrated with pipeline
    
    if (UNIPATH_BABY_STATS) {
// MakeDepend: dependency UnipathBabyStats
      make.AddRule( FilesIn( run_dir, "UnipathBabyStats.out" ),
		    JoinVecs( FilesIn( ref_dir, "genome.{fastb,lookup}",
				       "genome_extended.fasta"),
			      FilesIn( run_dir,"all_reads.{unibases,paths,paths_rc,"
				       "unipaths,pathsdb,unipathsdb}.kN",
				       "all_reads.unipaths.predicted_count.kN" ) ),
		    
		    "UnipathBabyStats " + pdrk
		    + ARG(NUM_THREADS, global_threads)
		    + ARG(READS, "all_reads"),
		    
		    "Evaluating unipaths for regional assembly" );
      
    }

    ///////////////////////////////////////////////////////////////////////////
    // Evaluate unipaths and copy number predictions
    //

// MakeDepend: dependency AddUnibaseJunctions
    if ( ADD_UNIBASE_JUNCTIONS ) {

      make.AddRule( FilesIn( run_dir, current_unipaths + ".junctionized.unibases.kN" ),
		    FilesIn( run_dir, current_unipaths + ".unibases.kN" ),
		    "AddUnibaseJunctions " + pdrk
		    + ARG(UNIBASES_IN, current_unipaths)
		    + ARG(UNIBASES_OUT, current_unipaths + ".junctionized"),
		    
		    "Add junctions to unibases." );
      
      current_unipaths += ".junctionized";

    }

    // Copy final unipaths.
    make.AddCpRule( run_dir + "/extended.fastb", run_dir + "/" +
		    current_unipaths + ".unibases.kN" );

    // Copy and create new filled reads kmerspace
    make.AddCpRule( run_dir, filledReads + "_ext",
		    run_dir, filledReads,
		    ".{fastb,pairs}");

      
    // Common path new unibases and filled reads
    AddRule_CommonPather(make, CP_THREADS, K, run_dir,
			 run_dir, "extended.fastb",
			 run_dir, filledReads + "_ext.fastb");

// MakeDepend: dependency MakeRcDb
    make.AddRule( FilesIn( run_dir, filledReads + "_ext.paths{_rc,db}.kN" ),
		  FilesIn( run_dir, filledReads + "_ext.paths.kN" ),
		  "MakeRcDb " + pdrk + 
		  " READS=" + filledReads + "_ext" );

// MakeDepend: dependency MakeRcDb
    make.AddRule( FilesIn( run_dir, "extended.paths{_rc,db}.kN" ),
		  FilesIn( run_dir, "extended.paths.kN" ),
		  "MakeRcDb " + pdrk + 
		  " READS=extended" );

// MakeDepend: dependency Unipather
    make.AddRule( FilesIn( run_dir, "extended.{unipaths{,db},unibases,unipath_adjgraph}.kN" ),
		  FilesIn( run_dir, "extended.{paths{,_rc,db}.kN,fastb}"),
		  "Unipather " + pdrk +
		  " READS=extended UNIPATHS=unipaths UNIBASES=unibases" +
		  " BUILD_UNIGRAPH=True",
		  
		  "Re-creating unipaths for extended unibases" );

// MakeDepend: dependency FilterPrunedReads
    if ( FILTER_PRUNED_READS )
    {
    make.AddRule( FilesIn( run_dir, filledReads + "_filt.{fastb,paths.kN,pairs,readtrack}" ),
		  FilesIn( run_dir, filledReads + "_ext.{fastb,paths.kN,pairs}",
			            "extended.unipathsdb.kN" ),
		  "FilterPrunedReads " + pdrk
		  + " READS_IN=" + filledReads + "_ext "
		  + " READS_OUT=" + filledReads + "_filt"
		  + " UNIPATHS=extended",
		  "Remove reads containing kmers removed during unipath extending/patching/etc." );
    }

    // update names to use later
    unibasesName = "extended";
    if ( FILTER_PRUNED_READS ) filledReads = filledReads + "_filt";
    else filledReads = filledReads + "_ext";
  } // endif USE_JUMPS_IN_UNIPATHS


  ///////////////////////////////////////////////////////////////////////////
  // Error correct jumping library reads
  //
  
  String jumpReadsInputName = "jump_reads_filt_cpd";

// MakeDepend: dependency CreateLookupTab
  make.AddRule( FilesIn( run_dir, unibasesName + ".unibases.kN.lookup" ),
		FilesIn( run_dir, unibasesName + ".unibases.kN"),
		"CreateLookupTab "
		" K=12"
		" IS_FASTB=True"
		" SOURCE=" + run_dir + "/" + unibasesName + ".unibases.kN"
		" OUTPUT=" + run_dir + "/" + unibasesName + ".unibases.kN.lookup",
		
		"Create lookup table from unibases" );


// MakeDepend: dependency ErrorCorrectJump
  make.AddRule( FilesIn( run_dir, "jump_reads_ec.{fastb,pairs,paths.kN,orig_id,readtrack}"),
		FilesIn( run_dir, jumpReadsInputName+ ".{fastb,qualb,pairs}",
			 unibasesName + ".{unibases,unipaths}.kN",
			 unibasesName + ".unibases.kN.lookup"),
		"ErrorCorrectJump  K=" + KS +
		//		  " PERFECT=True"
		" NUM_THREADS=" + ToString(ECJ_THREADS) +
		" QUERY_HEAD=" + run_dir + "/" + jumpReadsInputName +
		" MIN_MATCH=" + ToString(ECJ_MIN_MATCH) +
		" MISMATCH_THRESHOLD=" + ToString(ECJ_MISMATCH_THRESHOLD) +
		" MISMATCH_NEIGHBORHOOD=" + ToString(ECJ_MISMATCH_NEIGHBORHOOD) +
		" MISMATCH_BACKOFF=" + ToString(ECJ_MISMATCH_BACKOFF) +
		" SCORE_DELTA=" + ToString(ECJ_SCORE_DELTA) +
		" SCORE_MAX=" + ToString(ECJ_SCORE_MAX) +
		" OUT_HEAD=" + run_dir + "/jump_reads_ec "
		" REF_HEAD=" + run_dir + "/" + unibasesName,
		
		"Correct jumping library reads by aligning to fragment unibases." );


  if (USE_LONG_JUMPS) {

    if (ERROR_CORRECT_LONG_JUMPS) {
      
// MakeDepend: dependency ErrorCorrectJump
      make.AddRule( FilesIn( run_dir, "long_jump_reads_ec.{fastb,pairs,paths.kN,orig_id,readtrack}" ),
		    FilesIn( run_dir, long_jump_reads_current + ".{fastb,qualb,pairs}",
			     unibasesName + ".{unibases,unipaths}.kN",
			     unibasesName + ".unibases.kN.lookup"),
		    "ErrorCorrectJump  K=" + KS +
		    //		  " PERFECT=True"
		    " NUM_THREADS=" + ToString(ECJ_THREADS) +
		    " QUERY_HEAD=" + run_dir + "/" + "long_jump_reads_filt" +
		    " MIN_MATCH=" + ToString(ECJ_MIN_MATCH) +
		    " MISMATCH_THRESHOLD=" + ToString(ECJ_MISMATCH_THRESHOLD) +
		    " MISMATCH_NEIGHBORHOOD=" + ToString(ECJ_MISMATCH_NEIGHBORHOOD) +
		    " MISMATCH_BACKOFF=" + ToString(ECJ_MISMATCH_BACKOFF) +
		    " SCORE_DELTA=" + ToString(ECJ_SCORE_DELTA) +
		    " SCORE_MAX=" + ToString(ECJ_SCORE_MAX) +
		    " OUT_HEAD=" + run_dir + "/long_jump_reads_ec "
		    " REF_HEAD=" + run_dir + "/" + unibasesName +
		    " FLIP=False", // fosmid reads don't have to be reversed for correction
		    
		    "Correct long jumping library reads by aligning to fragment unibases." );
      
    } else {
      make.AddCpRule( run_dir, "long_jump_reads_ec", 
		      run_dir, long_jump_reads_current,
		      ".{fastb,pairs}" );
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Compute invariant size distributions of jump_reads_ec
	  

// MakeDepend: dependency SamplePairedReadDistributions
    make.AddRule( FilesIn( run_dir, "jump_reads_ec.distribs" ),
		  JoinVecs(  FilesIn( run_dir, "jump_reads_ec.{fastb,pairs}" ),
			     FilesIn( run_dir, unibasesName + ".unibases.kN" ) ),
		  String("SamplePairedReadDistributions") +
		  " NUM_THREADS=" + ToString(global_threads) +
		  " BUILD_LOOKUP_TABLE=True KEEP_LOOKUP_TABLE=False"
		  " UNIBASES_K=" + KS +
		  " UNIBASES=" + run_dir + "/" + unibasesName + ".unibases.k" + KS +
		  " HEAD_READS=" + run_dir + "/" + "jump_reads_ec"
		  " FLIP=False",
		  
		  "Use filled fragment unipaths to compute fragment read pair separations" );
 

  ///////////////////////////////////////////////////////////////////////////
  // Merge filled fragment and jumping reads together - no pathing required
  //

  String correctedJumpReads = "jump_reads_ec"; // not sure why this is needed twice

  if (MERGE_UNIPATHS_AS_READS) {

// MakeDepend: dependency SplitUnibases
    make.AddRule( FilesIn( run_dir, "unibases_as_reads.{fastb,pairs,"
			   "paths.kN}"),
		  FilesIn( run_dir, unibasesName + 
			   ".{unibases,unipaths,unipath_adjgraph}.kN"),
		  "SplitUnibases "
		  " IN_HEAD=" + run_dir + "/" + unibasesName +
		  " OUT_HEAD=" + run_dir + "/unibases_as_reads"
		  " K=" + KS,
		  
		  "Split the unibases into readlike objects" );

    AddRule_MergeReadSets(make, false, true,  true, false,    true,  false,
			  //    quals, pairs, paths, dist, tracker, repath //
			  K, CP_THREADS, run_dir,
			  "all_reads",  // target //
			  filledReads, "unibases_as_reads", correctedJumpReads);

  } else {

    AddRule_MergeReadSets(make, false, true,  true,  false,   true,  false,
			  //    quals, pairs, paths, dist, tracker, repath //
			  K, CP_THREADS, run_dir,
			  "all_reads",  // target //
			  filledReads, correctedJumpReads);

  }    
// MakeDepend: dependency MakeRcDb
  make.AddRule( FilesIn( run_dir, "all_reads.paths{_rc,db}.kN" ),
		FilesIn( run_dir, "all_reads.paths.kN" ),
		"MakeRcDb " + pdrk + 
		"READS=all_reads" );
  
  // Reuse extended unipaths rather than create new ones - they share the same kmerspace.
  make.AddCpRule( run_dir, "all_reads", 
		  run_dir, "extended",
		  ".{unipaths{,db},unibases,unipath_adjgraph}.kN" );
  
 
  ///////////////////////////////////////////////////////////////////////////
  // Find unipath copy number
  //

   
// MakeDepend: dependency UnibaseCopyNumber3
  make.AddRule( FilesIn( run_dir, "all_reads.unipaths.{predicted_count,kmer_hits,cn_raw}.kN" ),
		FilesIn( run_dir, "all_reads.unibases.kN", correctedFragReads + ".fastb",
			 correctedJumpReads + ".fastb", "ploidy"),
		
		"UnibaseCopyNumber3 " + pdrk +
		" NUM_THREADS=" + ToString(global_threads) +
		" READS=all_reads"
		" READS_EC=\"{" + correctedFragReads + "," + correctedJumpReads + "}\""
                + ARG(LOWER, UCN_LOWER) +
		" ERR_RATE=0.0 "
		" EXP_GAP_REMODEL=" + ( UCN_GAP_REMODEL ? "True" : "False" ),
		
		"Estimating the copy number of unibases (and the corresponding "
		"unipaths), based on how many k-mers pile on each unibase." );   


// MakeDepend: dependency UnipathEval
// Note that UnipathEval generates two PerfStats that are needed.
  if ( evalStandard ){
    make.AddOutputRule( FilesIn( run_dir, "all_reads.unipaths.kN.UnipathEval.out" ),
			JoinVecs( FilesIn( data_dir, "genome.{fastb,lookup}", "ploidy" ),
				  FilesIn( run_dir,"all_reads.{unipaths,unibases}.kN",
					   "all_reads.unipaths.predicted_count.kN") ),
			"UnipathEval " + pdrk +
			" READS=all_reads UNIPATHS=unipaths UNIBASES=unibases"
			" COPYNUM=True ARCHIVE=True "
			" NUM_THREADS=" + ToString(global_threads),
			
			"Evaluating how well we predicted the copy numbers of the unibases" );

    make.AddRule( FilesIn( run_dir, "all_reads.unipaths.kN.placements",
			   "all_reads.unipaths.k96.PathsToLocs.out"),
		  JoinVecs( FilesIn( data_dir, "genome.fastb" ),
			    FilesIn( run_dir,"all_reads.{paths,paths_rc,pathsdb}.kN",
				     "all_reads.{unipaths,unipathsdb,unibases}.kN",
				     "all_reads.unipaths.predicted_count.kN",
				     "all_reads.fastb") ),
		  "PathsToLocs " + pdrk +
		  " NUM_THREADS=" + ToString(CP_THREADS) +
		  " PATHS=all_reads.unipaths.kN" 
		  " READS=all_reads" 
		  " UNIPATHS_FASTB=all_reads.unibases.kN" 
		  " UNIPATH=True RENUMBER=True SHOW_PLACEMENTS=True"
		  " MATCH_LEN_PERCENT=100 RC_TOO=True ARCHIVE=True" 
		  " GENOME=genome"
		  " SAVE_PLACEMENTS_TO=" + run_dir + "/all_reads.unipaths.kN.placements",
		  
		  "Evaluating unipath coverage of a reference" );
  }


  ///////////////////////////////////////////////////////////////////////////
  // Align unibases to the reference genome to find their reference locations
  //
  
  if ( evalStandard ) {
// MakeDepend: dependency AlignUnibasesToRef
    make.AddRule( FilesIn( run_dir, "all_reads.unipaths.kN.locs" ),
		  JoinVecs( FilesIn( run_dir, "all_reads.unibases.kN" ),
			    FilesIn( data_dir, "genome.{fastb,lookup}" ) ),
		  "AlignUnibasesToRef " + pdrk +
		  " NUM_THREADS=" + ToString(global_threads) +
		  " READS=all_reads",
		  
		  "Finding locations of unibases on the reference" );
    
  }  // !evalStandard


  ///////////////////////////////////////////////////////////////////////////
  // Find possible read locations on the reference genome
  //

  if ( evalStandard ) {
    String readsAlignsFile, uniqueAlignsFile;
    AddRule_SQLAlign(make, "genome.fastb", "all_reads.fastb",
		     ref_dir, run_dir, "", 0, True, 0, 0,
		     "Aligning error corrected reads to the reference to generate"
		     "read locations");
    readsAlignsFile  = "all_reads.aligned.genome.sql.qltout";
    uniqueAlignsFile = "all_reads.aligned_unique.genome.sql.qltout";

// MakeDepend: dependency UniquifyAligns2
    make.AddRule( FilesIn( run_dir, uniqueAlignsFile ),
		  FilesIn( run_dir, "all_reads.pairs", readsAlignsFile),
		  "UniquifyAligns2 " + pdr + " READS=all_reads"
		  " ALL_ALIGNS_IN=" + readsAlignsFile + 
		  " UNIQ_ALIGNS_OUT=" + uniqueAlignsFile +
		  " MAX_DISCREP_TO_ACCEPT=3 MAX_DISCREP_TO_CONSIDER=4 "
		  " DATA_IN_RUN_DIR=True",
		  
		  "Finding better aligns of corrected reads to reference, based on read pairings" );

    if ( LR_READ_CLOUD_EVAL ) {
// MakeDepend: dependency GenomeReadLocsLG
      make.AddRule( FilesIn( run_dir, "all_reads.ref.locs" ),
		    JoinVecs(FilesIn( run_dir, "all_reads.fastb", uniqueAlignsFile),
			     FilesIn( data_dir, "genome.fastb") ),
		    
		    "GenomeReadLocsLG " + pdr +
		    " ALIGNS_IN=" + uniqueAlignsFile +
		    " READS=all_reads UNIQUE_ONLY=True",
		    
		    "Generates read location on genome file" );
    }
  }

  String current_scaffolds;
  String CLR_uni;

  if (BIG_MAP) {

// MakeDepend: dependency UnibaseCopyNumber3
  make.AddRule( FilesIn( run_dir, 
			 "extended.unipaths."
			 "{predicted_count,predicted_gap,kmer_hits,cn_raw}.kN" ),
		FilesIn( run_dir, "extended.unibases.kN", 
                         correctedFragReads + ".fastb", "ploidy"),
		
		"UnibaseCopyNumber3 " + pdrk +
		" NUM_THREADS=" + ToString(global_threads) +
		" READS=extended"
		" READS_EC=\"{" + correctedFragReads + "}\"" 
		" ERR_RATE=0.0 "
		" LOWER=True "
		" WRITE_GAP=True"
		" EXP_GAP_REMODEL=True",
		
		"Estimating the copy number of unibases (and the corresponding "
		"unipaths), based on how many k-mers pile on each unibase." );   

// MakeDepend: dependency CorrectLongReads
      make.AddRule( FilesIn( run_dir, "extended.long.right_exts" ),
		    FilesIn( run_dir, "extended.unibases.kN",
                             "frag_reads_edit.{fastb,qualb,pairs}",
			     "long_reads_orig.fastb"),

		    "CorrectLongReads " + pdrk
		    + ARG(NUM_THREADS, global_threads)
                    + ARG(KOUT, CLR_KOUT)
                    + ARG(USE_SHORTEST, True)
                    + ARG(NEW_FILTER, True)
                    + ARG(STANDARD_ALIGNS, True)
                    + ARG(FILTER, False)
                    + ARG(CLEAN_GRAPH, True)
		    + " IN_HEAD=extended",

		    "Correct long reads" );

// MakeDepend: dependency LongReadJoin
      make.AddRule( FilesIn( run_dir, "extended.long.{unipaths,unipathsdb,unibases,unipath_adjgraph}.k"
			     + ToString(CLR_KOUT) ),
		    FilesIn( run_dir, 
                         "{extended.long.right_exts,extended.unibases.kN{,.predicted_gaps.txt}}" ),

		    "LongReadJoin " + pdrk
		    + ARG(NUM_THREADS, global_threads)
                    + ARG(KBIG, CLR_KOUT)
		    + " IN_HEAD=extended",

		    "Join corrected long reads" );

// MakeDepend: dependency ShaveUnipathGraph
      make.AddRule( FilesIn( run_dir, "extended.long.shaved.{unipaths,unipathsdb,unibases,unipath_adjgraph}.k"
			     + ToString(CLR_KOUT)),
		    FilesIn( run_dir, "extended.long.{unipaths,unipathsdb,unibases,unipath_adjgraph}.k"
			     + ToString(CLR_KOUT)),
		  "ShaveUnipathGraph "
		  " DIR=" + run_dir +
		  " IN_HEAD=" + "extended.long"
		  " OUT_HEAD=" + "extended.long.shaved"
		  " K=" + ToString(CLR_KOUT) +
		  " NUM_THREADS=" + ToString(CP_THREADS),
		  
		  "Shave long unipaths");

    CLR_uni = "extended.long.shaved";

// MakeDepend: dependency FindUnipathGaps
    make.AddRule( FilesIn( run_dir, "extended.unibases.k" 
                       + ToString(K) + ".predicted_gaps.txt" ),
                  FilesIn( run_dir, "extended.unibases.kN",
                       "jump_reads_ec.{pairs,fastb}" ),
		  "FindUnipathGaps " + pdr +
                  ARG(K1, K) +
                  ARG(K, 40) +
		  " NUM_THREADS=" + ToString(global_threads) +
		  " HEAD=extended",
		  
		  "Compute gaps between unipaths" );

// MakeDepend: dependency FindUnipathGaps
    make.AddRule( FilesIn( run_dir, "extended.long.shaved.unibases.k640"
                       ".predicted_gaps.txt" ),
                  FilesIn( run_dir, "extended.long.shaved.unibases.k640",
                       "jump_reads_ec.{pairs,fastb}" ),
		  "FindUnipathGaps " + pdr +
                  ARG(K1, 640) +
                  ARG(K, 80) +
		  " NUM_THREADS=" + ToString(global_threads) +
		  " HEAD=extended.long.shaved",
		  
		  "Compute gaps between unipaths" );

// MakeDepend: dependency BigMap
    make.AddRule( FilesIn( sub_dir, "long_direct.{summary,superb{,.mapping},"
			   "contigs.{efasta,fasta,vecfasta,fastb,mapping},"
			   "assembly.{efasta,fasta}}" ),
		  FilesIn( run_dir, "extended.{unibases,unipaths{,.predicted_count,.cn_raw}}.kN",
			   "extended.long.shaved.unibases.k640",
                           "extended.unibases.k" + ToString(K) + ".predicted_gaps.txt",
                           "extended.long.shaved.unibases.k640.predicted_gaps.txt",
			   "jump_reads_ec.{fastb,pairs}"),
		  "BigMap " + pdr +
		  " NUM_THREADS=" + ToString(global_threads) +
		  " HEAD1=extended"
                  " VALIDATE=" + ( evalFull ? "True" : "False" )
		  + " HEAD2=extended.long.shaved",
		  
		  "Use BigMap to create genome map" );

    current_scaffolds = "long_direct";

  } else {


  ///////////////////////////////////////////////////////////////////////////
  // Build Unipath Link Graph
  //

// MakeDepend: dependency UnipathLocsLG
  make.AddRule( FilesIn( run_dir, "all_reads.unilocs." + KS + "{,.indexr}" ),
		FilesIn( run_dir, "all_reads.{paths{,_rc},unipaths{,db,.predicted_count}}.kN" ),
		"UnipathLocsLG " + pdrk +
		" NUM_THREADS=" + ToString(global_threads) +
		" READS=all_reads" + 
		" MAX_COPY_NUMBER=10 MIN_KMERS=1",
		  
		"For building a unipath link graph, aligning reads to unipaths" );


// MakeDepend: dependency BuildUnipathLinkGraphsLG
  make.AddRule( FilesIn( sub_dir, "unipath_link_graph.{seeds,cloud}.kN" ),
		FilesIn( run_dir, "{all_reads.{paths{,db,_rc},unipaths{,.predicted_count},"
			 "unipathsdb}.kN,all_reads.{unilocs." + KS + "{,.indexr},pairs}}",
			 "ploidy" ),
		"BuildUnipathLinkGraphsLG " + pdrsk +
		" READS=all_reads" +
		" MIN_UNIPATH_FOR_TRANSITIVE_JOIN=" 
		+ ToString(BULG_MIN_UNIPATH_FOR_TRANSITIVE_JOIN)
                + " NORMAL_MULT=" + ToString(BULG_LINK_GRAPH_NORMAL_MULT)
		+ " TRANSITIVE_FILL_IN=" 
		+ (BULG_TRANSITIVE_FILL_IN ? "True" : "False" ) +
		" NUM_THREADS=" + ToString(global_threads),
		
		"Building unipath link graph, showing approximate distance between pairs of unipaths" );

// MakeDepend: dependency EvalUnipathLinkGraphs
  if ( evalStandard && !BIG_MAP )
    make.AddRule( FilesIn( sub_dir, "EvalUnipathLinkGraphs.out" ),
		  JoinVecs( FilesIn( data_dir, "ploidy" ),
			    FilesIn( run_dir, "all_reads.{paths{,_rc}.kN,unipaths.kN,unipaths.kN.locs}" ),
			    FilesIn( sub_dir, "unipath_link_graph.seeds.kN" ) ),
		  "EvalUnipathLinkGraphs" + pdrsk
		  + " GRAPH=seeds READS=all_reads VERBOSE=True ",
		  
		  "Evaluate the unipath link graph by placing unipaths on reference" );


  ///////////////////////////////////////////////////////////////////////////
  // Select seed unipaths
  //

// MakeDepend: dependency SelectSeeds
  make.AddRule( FilesIn( sub_dir, "seeds.ids", "SelectSeeds.log"),
		JoinVecs( FilesIn( data_dir, (evalFull ? "genome.fastb" : "" ), "ploidy" ),
			  FilesIn( run_dir, "all_reads.{pairs,{paths,paths_rc,pathsdb}.kN,"
				   "unipaths{,db,.predicted_count}.kN,"
				   "{unilocs." + KS + ",pairs}}",
				   ( evalFull ? "all_reads.unipaths.kN.locs" : "") ),
			  FilesIn( sub_dir, "unipath_link_graph.{cloud}.kN" ) ),
		"SelectSeeds " + pdrsk +
		" NUM_THREADS=" + ToString(global_threads) +
		" READS=all_reads"
                " CHEAT_RC=" + ToStringBool(EVALUATION == "CHEAT") +
		" EVAL=" + ToStringBool(evalFull),
		
		"Select seed unipaths" );

  if (evalStandard && !BIG_MAP) {
// MakeDepend: dependency MapSeeds
    make.AddRule( FilesIn( sub_dir, "seeds.map" ),
		  JoinVecs( FilesIn( sub_dir, "{seeds.ids,unipath_link_graph.seeds.kN}"),
			    FilesIn( data_dir, "genome.lookup"),
			    FilesIn( run_dir, "all_reads.{unibases,unipaths}.kN") ),
		  "MapSeeds " + pdrsk +
		  " NUM_THREADS=" + ToString(global_threads) +
		  " READS=all_reads",
		  
		  "Find true location of seeds on the reference - creates seeds.map" );
  }


  ///////////////////////////////////////////////////////////////////////////
  // Run LocalizeReadsLG
  //

  vec<filenamepart_t> LocalizeReadsDeps = 
    JoinVecs( FilesIn( data_dir, (evalFull ? "genome.fastb" : "" ), "ploidy" ),
	      FilesIn( run_dir, "all_reads.{fastb,pairs}",
		       "all_reads.{paths{,_rc,db},unibases,unipaths{,db,.predicted_count}}.kN",
		       "all_reads.unilocs." + KS,
		       ( evalFull ? "all_reads.unipaths.kN.locs" : ""),
		       (LR_READ_CLOUD_EVAL ? "all_reads.ref.locs" : "") ),
	      FilesIn( sub_dir, "unipath_link_graph.cloud.kN", "seeds.ids" ) );
  
// MakeDepend: dependency LocalizeReadsLG
  make.AddRule( FilesIn( sub_dir, "LocalizeReadsLG.log", "localized.hbvs" ),
		LocalizeReadsDeps,
		"LocalizeReadsLG " + pdrsk +
		" READS=all_reads" +
		" NUM_THREADS=" + ToString(LR_THREADS) +
		(evalFull ? " USE_TRUTH=True" : " USE_TRUTH=False") +
		" EVAL=" + ToStringBool(LR_EVAL) +
		" READ_CLOUD_EVAL=" + ToStringBool( LR_READ_CLOUD_EVAL ) +
		ARG(LOCAL_DUMP, LR_LOCAL_DUMP) +
		ARG(OVERRIDE_TIMEOUTS, LR_OVERRIDE_TIMEOUTS) +
		ARG(BRIDGE_BETWEEN_CLOUD_UNIPATHS,LR_BRIDGE_BETWEEN_CLOUD_UNIPATHS) +
		" NHOOD_GLUEPERFECT_SIZE=" + ToString( LR_NHOOD_GLUEPERFECT_SIZE ) +
		" MAX_COPY_NUMBER_OTHER=" + ToString( LR_NHOOD_CN_FACTOR * ploidy ) +
		" " + LR_PASS_THRU,
		
		"Performing an ALLPATHS assembly from the read-estimated unipaths "
		"and the error-corrected paired reads" );
  
  
  ///////////////////////////////////////////////////////////////////////////
  // Merge neighborhoods
  //

  make.AddCpRule( sub_dir, "reads", run_dir, "all_reads", ".unipath_adjgraph.kN" );
  
// MakeDepend: dependency MergeNeighborhoods1
  make.AddRule( FilesIn( sub_dir, "reads.{fastb,{paths,paths_rc,pathsdb}.kN}",
			 "reads.{unipaths,unibases}.kN", "nhood.hypers" ),
		JoinVecs( FilesIn( run_dir, "all_reads.{unibases}.kN" ),
			  FilesIn( sub_dir, "seeds.ids", "localized.hbvs" ) ),
		"MergeNeighborhoods1 " + pdrsk
		+ " NUM_THREADS=" + ToString( CP_THREADS ),
		"Performing an ALLPATHS assembly from the read-estimated unipaths "
		"and the error-corrected paired reads, phase 1" );
  
// MakeDepend: dependency MergeNeighborhoods2
  make.AddRule( FilesIn( sub_dir, "hyper.prelim{,.dot}" ),
                JoinVecs( FilesIn( data_dir, "ploidy" ),
			  FilesIn( run_dir, "all_reads.{unibases,unipaths.predicted_count}.kN" ),
			  FilesIn( sub_dir, "reads.paths{,db}.kN", "nhood.hypers" ) ),
		"MergeNeighborhoods2 " + pdrsk
		+ ARG(NUM_THREADS, MN_THREADS) 
                + ARG(MIN_OVERLAP, MN_MIN_OVERLAP)
                + ARG(MIN_PROPER_OVERLAP, MN_MIN_PROPER_OVERLAP),
		
		"Performing an ALLPATHS assembly from the read-estimated unipaths "
		"and the error-corrected paired reads, phase 2" );

// MakeDepend: dependency MergeNeighborhoods3
  make.AddRule( FilesIn( sub_dir, "hyper{,.dot}" ),
                FilesIn( sub_dir, "hyper.prelim" ),
		"MergeNeighborhoods3 " + pdrsk,
		"Performing an ALLPATHS assembly from the read-estimated unipaths "
		"and the error-corrected paired reads, phase 3" );


  ///////////////////////////////////////////////////////////////////////////
  // Recover Missing Unipaths
  //

// MakeDepend: dependency RecoverUnipaths
  make.AddRule( JoinVecs( FilesIn( sub_dir, "hyper.extra.fastb",
				   "hyper.recovered{,.fasta,.fastb,.graphml}",
				   (evalFull ? "hyper.extra.to_align" : "") ),
			  FilesIn( sub_dir + "/recover", "reads.{fastb,paths{,_rc,db}.kN}" ) ),
		      JoinVecs(FilesIn( data_dir, (evalFull ? "genome.lookup" : "") ),
			       FilesIn( sub_dir, "reads.{paths{,_rc,db},unipaths,unibases}.kN",
					"reads.{unipath_adjgraph.kN,fastb}", "hyper" ) ),
		      "RecoverUnipaths " + pdrsk 
		      + ARG(NUM_THREADS, global_threads) 
                      + ARG(MIN_KEEPER, MIN_CONTIG) +
		      " WRITE=True "
		      " READS=reads"
		      " WRUN_IN="
		      " WRUN_OUT=recover"
		      " HYPER_OUT=hyper.recovered"
		      " USE_TRUTH=" + ToStringBool(evalFull),
		      
		      "recover unipaths that did not make it into the assembly." );


  make.AddCpRule( sub_dir + "/hyper_plus", sub_dir + "/hyper.recovered" );

  ///////////////////////////////////////////////////////////////////////////
  // Write graph edges as a fasta
  //

// MakeDepend: dependency DumpHyper
  make.AddRule( FilesIn( sub_dir, "hyper.{fasta,fastb,graphml}" ),
		FilesIn( sub_dir, "hyper", "reads.{fastb,paths{,_rc,db}.kN}" ),
		"DumpHyper " + pdr +
		" SUBDIR=" + SUBDIR + " WRUN=" );
  

  ///////////////////////////////////////////////////////////////////////////
  // Create Scaffold
  //

  // Gather together reads for scaffolding

  if (USE_LONG_JUMPS)
// MakeDepend: dependency MergeReadSets
    make.AddRule( FilesIn( run_dir, "scaffold_reads.{fastb,pairs,distribs,readtrack}"),
		  FilesIn( run_dir, "long_jump_reads_ec.{fastb,pairs,readtrack}",
			   correctedJumpReads + ".{fastb,pairs,distribs,readtrack}"),
		  "MergeReadSets "
		  " DIR=" + run_dir +
		  " READS_IN=" + "\"{"+ correctedJumpReads + ",long_jump_reads_ec}\""
		  " READS_OUT=scaffold_reads"
		  " K=" + ToString(K) +
		  " NUM_THREADS=" + ToString(CP_THREADS) +
		  " MERGE_QUALS=False MERGE_PATHS=False REPATH=False"
		  " MERGE_PAIRS=True TRACK_READS=True"
		  " MERGE_DISTRIBS=True FORCE_MERGE_DISTRIBS=True",
		
		"Merging scaffolding reads");

  else {
    make.AddCpRule( run_dir, "scaffold_reads", 
		    run_dir, correctedJumpReads,
		    ".{fastb,pairs,distribs}");
  }


// MakeDepend: dependency FlattenHKP
  String FlattenHKP_target = ASSISTED_PATCHING ? "initial_scaffolds0" : "initial_scaffolds";
  make.AddRule( FilesIn( sub_dir, "hyper_plus.FlattenHKP.log", FlattenHKP_target + 
			       ".{superb,contigs.fasta,assembly.fasta}"),
		      JoinVecs( FilesIn( sub_dir, "hyper_plus" ),
				FilesIn( sub_dir + "/recover", 
					 "reads.{fastb,paths{,_rc,db}.kN}"  )),
		      "FlattenHKP " + pdrs + 
		      " HYPER=hyper_plus"
		      + ARG(SCAFFOLDS, FlattenHKP_target) 
		      + ARG(NUM_THREADS, global_threads) 
                      + ARG(MIN_EDGE_TO_SAVE, MIN_CONTIG)
                      + ARG(NEW_ALGORITHM, FLA_NEW_ALGORITHM) + " WRUN=recover",
		      
		      "Flatten the assembly HKP into fastavectors, to make a trivial scaffold structure." );
  
  if (ASSISTED_PATCHING) {
    String current_scaffolds = FlattenHKP_target;
    String assembly_patches;
    
    // MakeDepend: dependency ScaffoldContigsOnReference
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".ref.{superb,contigs.fasta,contigs.fastb}"),
		  JoinVecs(FilesIn( sub_dir,  current_scaffolds + ".{superb,contigs.fasta}"),
			   FilesIn( ref_dir,  "genome.lookup")),
		  "ScaffoldContigsOnReference "
		  " HEAD_REF=" + ref_dir + "/genome"
		  " ASSEMBLY_IN=" + sub_dir + "/" + current_scaffolds +
		  " ASSEMBLY_OUT=" + sub_dir + "/" + current_scaffolds + ".ref",
		  "Use reference for initial contig layout for patching" );

    current_scaffolds += ".ref";

    // MakeDepend: dependency FastaToEfasta
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".contigs.efasta"),
		  FilesIn( sub_dir, current_scaffolds + ".contigs.fasta"),
		  "FastaToEfasta "
		  " IN=" + sub_dir + "/" + current_scaffolds + ".contigs.fasta"
		  " OUT=" + sub_dir + "/" + current_scaffolds + ".contigs.efasta",
		  "Convert reference-aligned initial contigs to efasta" );
    
    if (PATCH_SCAFFOLDS) {
      // MakeDepend: dependency PostPatcher
      make.AddRule( JoinVecs( FilesIn( sub_dir, current_scaffolds + ".patched.edits" ),
			      FilesIn( sub_dir, "post_patch/PostPatcher." + current_scaffolds + "."
				       "{INGAP,JDLEN,JINGAP,JOINDATA,TIGS,patched.RESULTS}" )),
		    JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				       "{superb,contigs.fasta}" ),
			      FilesIn( run_dir, final_uni_for_UnipathPatcher + ".unibases.kN",
				       "frag_reads_filt_cpd.{fastb,qualb,pairs}",
				       "jump_reads_filt_cpd.{fastb,qualb,pairs}" ),
			      FilesIn( run_dir + "/unipath_patch", "UnipathPatcher.{SEGS,JSEGS}") ),
		    "PostPatcher " + pdrsk 
		    + ARG(NUM_THREADS, global_threads)
		    + ARG(SCAFFOLDS_IN, current_scaffolds)
		    + ARG(FRAG_READS, "frag_reads_filt_cpd")
		    + ARG(JUMP_READS, "jump_reads_filt_cpd")
		    + ARG(MAX_MEMORY_GB, MAX_MEMORY_GB)
		    + ARG(MOC_STRINGENT, PP_MOC_STRINGENT)
		    + ARG(HEAD, final_uni_for_UnipathPatcher),
		    
		    "Patch gaps in initial scaffolds." );
      
      assembly_patches += (assembly_patches != "" ? "," : "") 
	+ current_scaffolds + ".patched.edits";
    }
      
    
    if (LONG_READ_POST_PATCH) {
      // MakeDepend: dependency LongReadPostPatcher
      make.AddRule( FilesIn( sub_dir, current_scaffolds + ".longread.edits" ),
		    JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				       "{superb,contigs.{fasta,fastb}}" ),
			      FilesIn( run_dir, "long_reads_orig.fastb" ) ),
		    "LongReadPostPatcher " + pdrs
		    + ARG(NUM_THREADS, global_threads)
		    + ARG(SCAFFOLDS_IN, current_scaffolds)
		    + ARG(SCAFFOLDS_OUT, current_scaffolds + ".longread")
		    + ARG(READS, "long_reads_orig"),
		  
		    "Patch gaps in initial scaffolds using long reads." );
      
      // MakeDepend: dependency CleanLongReadPatches
      make.AddRule( FilesIn( sub_dir, current_scaffolds + ".longread.fixed.edits" ),
		    JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				       "{longread.edits,superb,contigs.fasta}" ),
			      FilesIn( run_dir, "all_reads.unibases.kN" ) ),
		    "CleanLongReadPatches " + pdrsk
		    + ARG( NUM_THREADS, global_threads) 
		    + ARG(SCAFFOLDS_IN, current_scaffolds)
		    + ARG(EDITS_IN, current_scaffolds + ".longread.edits")
		    + ARG(EDITS_OUT, current_scaffolds + ".longread.fixed.edits"),
		  
		    "Clean up long read patches in initial scaffolds." );
      
      assembly_patches += (assembly_patches != "" ? "," : "") 
	+ current_scaffolds + ".longread.fixed.edits";
      
    }


    if ( PATCH_SCAFFOLDS || LONG_READ_POST_PATCH ) {
      // MakeDepend: dependency ApplyGapPatches
      make.AddRule( FilesIn( sub_dir, current_scaffolds + ".applied."
			     "{contigs.{fasta,fastb,efasta,mapping,vecfasta},"
			     "assembly.{fasta,efasta},superb{,.mapping},summary}" ),
		    JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				       "{superb,contigs.{fasta,efasta}}" ),
			      FilesIn( sub_dir, "{" + assembly_patches + "}" ) ),
		    "ApplyGapPatches " + pdrs
		    + ARG(SCAFFOLDS_IN, current_scaffolds)
		    + ARG( EDITS_IN, "\"{" + assembly_patches + "}\"")
		    + ARG(SCAFFOLDS_OUT, current_scaffolds + ".applied"),
		    
		    "Patch gaps in initial scaffolds." );
      
      current_scaffolds += ".applied";
    }

    make.AddCpRule( sub_dir + "/initial_scaffolds.contigs.fasta",
		    sub_dir + "/" + current_scaffolds + ".contigs.fasta" );
    
    // MakeDepend: dependency TrivialSupersFromContigs
    make.AddRule( FilesIn( sub_dir, "initial_scaffolds.superb"),
		  FilesIn( sub_dir, "initial_scaffolds.contigs.fasta"),
		  "TrivialSupersFromContigs "
		  " HEAD=" + sub_dir + "/initial_scaffolds",
		  
		  "Re-generate trivial initial scaffolds (singletons from contigs)" );
  } // end assisted patching block


// MakeDepend: dependency AlignPairsToFasta
  make.AddRule( FilesIn( sub_dir, "scaffold_reads.qltoutlet{,.index}" ),
                JoinVecs( FilesIn( sub_dir, "initial_scaffolds.contigs.fasta" ),
                          FilesIn( run_dir, "scaffold_reads.{fastb,pairs}" ) ),
                "AlignPairsToFasta " + pdrs +
		" NUM_THREADS=" + ToString( global_threads ) +
		" WRUN=recover READS=scaffold_reads ALIGNS=scaffold_reads"
		" FASTA=initial_scaffolds.contigs",

		"Align jumping reads to the assembly fasta." );

// MakeDepend: dependency RemoveHighCNAligns
  make.AddRule( FilesIn( sub_dir, "scaffold_reads_filtered.qltoutlet.index" ),
		JoinVecs( FilesIn( sub_dir, "scaffold_reads.qltoutlet{,.index}",
				   "initial_scaffolds.contigs.fasta" ),
			  FilesIn( run_dir, "all_reads.{unibases,unipaths.predicted_count}.k" + KS,
				   "ploidy") ),
		"RemoveHighCNAligns " + pdrsk +
		" NUM_THREADS=" + ToString( CP_THREADS ) +
		" UNIPATHS=all_reads"
		" ALIGNS_IN=scaffold_reads"
		" ALIGNS_OUT=scaffold_reads_filtered"
		" SCAFFOLDS_IN=initial_scaffolds"
		" FORCE=True",
		
		"Remove reads that align portion of contigs with CN greater than ploidy" );
  
  make.AddCpRule( sub_dir + "/scaffold_reads_filtered.qltoutlet",
		  sub_dir + "/scaffold_reads.qltoutlet");
  
// MakeDepend: dependency MakeScaffoldsLG
  make.AddRule( FilesIn( sub_dir, "linear_scaffolds0.{superb,contigs.fasta,assembly.fasta,summary,is_circular}" ),
		JoinVecs( FilesIn( sub_dir, "{scaffold_reads_filtered.qltoutlet{,.index},initial_scaffolds.{superb,contigs.fasta}}" ),
			  FilesIn( run_dir, "scaffold_reads.{fastb,pairs}", "ploidy" ),
			  FilesIn( run_dir, MS_FIT_GAPS? "scaffold_reads.distribs":"" ),
			  FilesIn( run_dir, MS_USE_UNIBASES? "all_reads.{unibases,unipaths.predicted_count}.kN":"") ),
		"MakeScaffoldsLG " + pdrs +
		" REGAP_NEG_GAPS=" + ToStringBool( MS_REGAP_NEG_GAPS ) +
		" REGAP=" + ToStringBool( MS_REGAP ) +
		" READS=scaffold_reads ALIGNS=scaffold_reads_filtered"
		" SCAFFOLDS_IN=initial_scaffolds HEAD_OUT=linear_scaffolds0 "
		+ ARG(NUM_THREADS, global_threads) 
                + ARG(FIT_GAPS, MS_FIT_GAPS) 
                + ARG(USE_UNIBASES, MS_USE_UNIBASES)
		+ ARG(UNIBASES_K, K)
		+ ARG(UNIBASES_HEAD, "all_reads")
                // + ARG(FIT_GAPS_EXP, MS_FIT_GAPS_EXP) 
                + " " + MS_PASS_THRU,
		
		"Use linking info to build true scaffolds out of a trivial scaffold structure." );

  current_scaffolds = "linear_scaffolds0";


} // BIG MAP V LRLG 

  if (CLEAN_ASSEMBLY) {
// MakeDepend: dependency CleanAssembly
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".clean.{contigs.{efasta,fasta,fastb,vecfasta,mapping},"
			   "assembly.{fasta,efasta},superb{,.mapping},summary}"),
		  FilesIn( sub_dir, current_scaffolds + ".{contigs.fasta,superb}" ),
		  "CleanAssembly "
		  " HEAD_IN=" + sub_dir + "/" + current_scaffolds +
		  " HEAD_OUT=" + sub_dir + "/" + current_scaffolds + ".clean"
		  " MIN_SCAFFOLD_SIZE=1000"
		  " MIN_CONTIG_SIZE_SOLO=0"
		  + ARG(NUM_THREADS, global_threads)
                  + ARG(MIN_CONTIG_SIZE_IN, CA_MIN_CONTIG_SIZE_IN)
                  + ARG(MIN_INDIRECT_GAP, ( LONG_READ_POST_PATCH ?
                          CA_MIN_INDIRECT_GAP_FOR_HYBRID : 0 ))
                  + ARG(MIN_UNIQUE, CA_MIN_UNIQUE)
                  + ARG(MAX_SWALLOW, CA_MAX_SWALLOW),
		  
		  "Remove small scaffolds (and/or contigs) from assembly" );
    
    current_scaffolds += ".clean";
  } else {
    make.AddCpRule( sub_dir + "/" + current_scaffolds + ".contigs.efasta",
		    sub_dir + "/" + current_scaffolds + ".contigs.fasta");
  }

  if (REMODEL) {
// MakeDepend: dependency RemodelGaps
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".remodel." +
			   "{contigs.{efasta,fasta,fastb,vecfasta},"
			   "superb,assembly.{fasta,efasta},summary}" ),
		  JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				     "{contigs.efasta,contigs.fastb,superb}" ),
			    FilesIn( run_dir,
                                     "jump_reads_ec.{fastb,pairs}" ) ),

		  "RemodelGaps " + pdrs + ARG(SCAFFOLDS_IN, current_scaffolds)
		  + ARG(SCAFFOLDS_OUT, current_scaffolds + ".remodel")
		  + ARG( NUM_THREADS, global_threads ) 
		  + ARG(JUMP_READS_IN, "jump_reads_ec")
                  + ARG(WRITE, True),

		  "Recompute some gaps.");

    current_scaffolds += ".remodel";
  }

  String assembly_patches;

  if (PATCH_SCAFFOLDS) {
// MakeDepend: dependency PostPatcher
    make.AddRule( JoinVecs( FilesIn( sub_dir, current_scaffolds + ".patched.edits" ),
			    FilesIn( sub_dir, "post_patch/PostPatcher." + current_scaffolds + "."
				     "{INGAP,JDLEN,JINGAP,JOINDATA,TIGS,patched.RESULTS}" )),
		  JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				     "{superb,contigs.fasta}" ),
			    FilesIn( run_dir, final_uni_for_UnipathPatcher + ".unibases.kN",
				     "frag_reads_filt_cpd.{fastb,qualb,pairs}",
				     "jump_reads_filt_cpd.{fastb,qualb,pairs}" ),
			    FilesIn( run_dir + "/unipath_patch", "UnipathPatcher.{SEGS,JSEGS}") ),
		  "PostPatcher " + pdrsk 
		  + ARG(NUM_THREADS, global_threads)
		  + ARG(SCAFFOLDS_IN, current_scaffolds)
		  + ARG(FRAG_READS, "frag_reads_filt_cpd")
		  + ARG(JUMP_READS, "jump_reads_filt_cpd")
		  + ARG(MAX_MEMORY_GB, MAX_MEMORY_GB)
                  + ARG(MOC_STRINGENT, PP_MOC_STRINGENT)
		  + ARG(HEAD, final_uni_for_UnipathPatcher),
		  
		  "Patch gaps in scaffolds." );

  assembly_patches += (assembly_patches != "" ? "," : "") 
       + current_scaffolds + ".patched.edits";

  }

  if (LONG_READ_POST_PATCH) {
// MakeDepend: dependency LongReadPostPatcher
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".longread.edits" ),
		  JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				     "{superb,contigs.{fasta,fastb}}" ),
			    FilesIn( run_dir, "long_reads_orig.fastb" ) ),
		  "LongReadPostPatcher " + pdrs
		  + ARG(NUM_THREADS, global_threads)
		  + ARG(SCAFFOLDS_IN, current_scaffolds)
		  + ARG(SCAFFOLDS_OUT, current_scaffolds + ".longread")
		  + ARG(READS, "long_reads_orig"),

		  
		  "Patch gaps in scaffolds using long reads." );

// MakeDepend: dependency CleanLongReadPatches
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".longread.fixed.edits" ),
		  JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				     "{longread.edits,superb,contigs.fasta}" ),
			    FilesIn( run_dir, "all_reads.unibases.kN" ) ),
		  "CleanLongReadPatches " + pdrsk
		  + ARG( NUM_THREADS, global_threads) 
                  + ARG(SCAFFOLDS_IN, current_scaffolds)
		  + ARG(EDITS_IN, current_scaffolds + ".longread.edits")
		  + ARG(EDITS_OUT, current_scaffolds + ".longread.fixed.edits"),
		  
		  "Clean up long read patches." );

  assembly_patches += (assembly_patches != "" ? "," : "") 
       + current_scaffolds + ".longread.fixed.edits";

  }

  if ( KPATCH && BIG_MAP ) {
// MakeDepend: dependency KPatch
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".kpatch.edits" ),
		  JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				     "{contigs.{efasta,fastb},superb}" ),
			    FilesIn( run_dir,  
				     ( !BIG_MAP ? "all_reads.unibases.k" + KS
                                 : CLR_uni + ".unibases.k" + ToString(CLR_KOUT)
                                      ) ) ),
		  "KPatch " + pdrs
                  + ARG(K, ( !BIG_MAP ? K : CLR_KOUT ))
                  + ARG(HEAD, ( !BIG_MAP ? "all_reads" : CLR_uni ))
		  + ARG(JUMPS_IN, "jump_reads_ec")
		  + ARG(SCAFFOLDS_IN, current_scaffolds)
		  + ARG(SCAFFOLDS_OUT, current_scaffolds + ".kpatch"),

		  "Patch some gaps in the final assembly by walking through "
                  "the unipath graph");

      assembly_patches += (assembly_patches != "" ? "," : "") 
	+ current_scaffolds + ".kpatch.edits";

  }


  if ( PATCH_SCAFFOLDS || LONG_READ_POST_PATCH ) {
// MakeDepend: dependency ApplyGapPatches
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".applied."
			   "{contigs.{fasta,fastb,efasta,mapping,vecfasta},"
			   "assembly.{fasta,efasta},superb{,.mapping},summary}" ),
		  JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				     "{superb,contigs.{fasta,efasta}}" ),
			    FilesIn( sub_dir, "{" + assembly_patches + "}" ) ),
		  "ApplyGapPatches " + pdrs
		  + ARG(SCAFFOLDS_IN, current_scaffolds)
		  + ARG( EDITS_IN, "\"{" + assembly_patches + "}\"")
		  + ARG(SCAFFOLDS_OUT, current_scaffolds + ".applied"),
		  
		  "Apply gap patches." );

   current_scaffolds += ".applied";

   }

  if (CONNECT_SCAFFOLDS) {
// MakeDepend: dependency AlignReads
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".readlocs" ),
		  JoinVecs( FilesIn( run_dir, "frag_reads_filt_cpd.{fastb,pairs}",
                                              "jump_reads_filt_cpd.{fastb,pairs}"
                                              /*
                                                  ,
                                              "long_jump_reads_filt.{fastb,pairs}" 
                                              */
                                                  ),
			    FilesIn( sub_dir, 
                                 current_scaffolds + ".{contigs.fastb,superb}" ) ),
                  "AlignReads " + pdrs 
		  + ARG( NUM_THREADS, global_threads) 
		  + ARG( ASSEMBLY, current_scaffolds ),

                  "Build read locations on the assembly for some of the reads." );

// MakeDepend: dependency ConnectScaffolds
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".connected.{contigs.fasta,contigs.vecfasta,contigs.fastb,assembly.fasta,superb}" ),
		  FilesIn( sub_dir, current_scaffolds + ".{readlocs,contigs.fasta,contigs.fastb,superb}" ),
                  "ConnectScaffolds " + pdrs + ARG( ASSEMBLY, current_scaffolds ),

                  "Connect some scaffolds that overlap." );

   current_scaffolds += ".connected";
  }

// MakeDepend: dependency AlignReads
  make.AddRule( FilesIn( sub_dir, current_scaffolds + ".readlocs" ),
		  JoinVecs( FilesIn( run_dir, "frag_reads_filt_cpd.{fastb,pairs}",
                                              "jump_reads_filt_cpd.{fastb,pairs}"
                                              /*
                                                  ,
                                              "long_jump_reads_filt.{fastb,pairs}" 
                                              */
                                               ),
			  FilesIn( sub_dir, 
				   current_scaffolds + ".{contigs.fastb,superb}" ) ),
		"AlignReads " + pdrs 
		+ ARG( NUM_THREADS, global_threads) 
		+ ARG( ASSEMBLY, current_scaffolds ),
		
		"Build read locations on the assembly for some of the reads." );

  string aligned_scaffolds = current_scaffolds;

// MakeDepend: dependency GenerateMapUnibasesToContigs
  make.AddRule( FilesIn( sub_dir, current_scaffolds + ".u2c" ),
		JoinVecs( FilesIn( run_dir, "all_reads.unibases.k" + KS ),
			  FilesIn( sub_dir, current_scaffolds + ".contigs.fastb" ) ),
		"GenerateMapUnibasesToContigs " + pdrsk
		+ ARG( NUM_THREADS, global_threads) 
		+ ARG( ASSEMBLY, current_scaffolds ),

		"Generate a unibases to contigs map (saved as a UInt64VecVec)" );

// MakeDepend: dependency TagCircularScaffolds
  make.AddRule( FilesIn( sub_dir, current_scaffolds + ".rings" ),
		JoinVecs( FilesIn( run_dir, "frag_reads_filt_cpd.pairs",
				   "jump_reads_filt_cpd.pairs",
				   (USE_LONG_JUMPS ? "long_jump_reads_filt.pairs" : "")),
			  FilesIn( sub_dir, current_scaffolds + ".{contigs.{fasta},superb}",
				   aligned_scaffolds + ".readlocs") ),
		"TagCircularScaffolds " + pdrs 
		+ ARG( NUM_THREADS, global_threads) +
		" FRAG_READS=frag_reads_filt_cpd JUMP_READS=jump_reads_filt_cpd"
		+ ARG( LONG_JUMP_READS, (USE_LONG_JUMPS ? "long_jump_reads_filt" : ""))
		+ ARG( ASSEMBLY, current_scaffolds ),
		
		"Identify and tag circular scaffolds" );
  
  bool have_efasta = false;

  String assembly_edits;

  if ( FIX_SOME_INDELS || FIX_ASSEMBLY_BASE_ERRORS ) {
// MakeDepend: dependency FixPrecompute
    make.AddRule( FilesIn( sub_dir, 
                       "post_patch/PostPatcher." + current_scaffolds + "."
                       + "{TIGS,JRALIGNS,UALIGNS,JRALIGNS_PSORT,JRPS_START}",
                       "post_patch/PostPatcher.stable.JREADLEN" ),
		  JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				     "{contigs.fasta,contigs.vecfasta}" ),
			    FilesIn( run_dir, "extended40.shaved.unibases.k" + KS,
				     "jump_reads_filt_cpd.fastb" ),
			    FilesIn( run_dir + "/unipath_patch", "UnipathPatcher.JSEGS") ),
		  "FixPrecompute " + pdrsk
		  + ARG( NUM_THREADS, global_threads) 
		  + ARG( SCAFFOLDS_IN, current_scaffolds)
		  + ARG( JUMP_READS, "jump_reads_filt_cpd")
		  + ARG( HEAD, final_uni_for_UnipathPatcher),

		  "Precompute stuff for FixSomeIndels and FixAssemblyBaseErrors");

  }

  if (FIX_SOME_INDELS) {
// MakeDepend: dependency FixSomeIndels
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".indel.edits" ),
		  JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				     "{contigs.vecfasta,readlocs}",
				     "post_patch/PostPatcher." + current_scaffolds + "."
				     + "{JRALIGNS_PSORT,JRPS_START}"),
			    FilesIn( run_dir,"frag_reads_filt.{fastb,qualb}",
				     "jump_reads_filt_cpd.{fastb,qualb}" ) ),
		  "FixSomeIndels " + pdrsk
		  + ARG( NUM_THREADS, global_threads) +
		  " WRITE=False " 
		  + ARG( SCAFFOLDS_IN, current_scaffolds)
		  + ARG( JUMP_READS, "jump_reads_filt_cpd")
		  + ARG( HEAD, final_uni_for_UnipathPatcher)
		  + ARG( SCAFFOLDS_OUT, current_scaffolds + ".indel"),

		  "Fix minor indel errors in final assembly");


    assembly_edits += (assembly_edits != "" ? "," : "") +  current_scaffolds + ".indel.edits";
  }

  if (FIX_ASSEMBLY_BASE_ERRORS) {
// MakeDepend: dependency FixAssemblyBaseErrors
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".base.edits"),
		  JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				     "{contigs.fasta,contigs.vecfasta,superb,readlocs}",

                       "post_patch/PostPatcher." + current_scaffolds + "."
                       + "{TIGS,JRALIGNS,UALIGNS,JRALIGNS_PSORT,JRPS_START}",
                       "post_patch/PostPatcher.stable.JREADLEN" ),

			    FilesIn( run_dir, "extended40.shaved.unibases.k" + KS,
				     "frag_reads_filt_cpd.{fastb,qualb,pairs}",
				     "jump_reads_filt_cpd.{fastb,qualb,pairs}" ) ),
		  "FixAssemblyBaseErrors " + pdrsk +
		  " WRITE=False " 
		  + ARG( NUM_THREADS, global_threads)
		  + ARG( SCAFFOLDS_IN, current_scaffolds)
		  + ARG( JUMP_READS, "jump_reads_filt_cpd")
		  + ARG( HEAD, final_uni_for_UnipathPatcher)
		  + ARG( SCAFFOLDS_OUT, current_scaffolds + ".base"),

		  "Fix minor base errors in final assembly");

    assembly_edits += (assembly_edits != "" ? "," : "") +  current_scaffolds + ".base.edits";
  }

  if (FIX_SOME_INDELS || FIX_ASSEMBLY_BASE_ERRORS) {
// MakeDepend: dependency ApplyAssemblyEdits
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".fixed." +
			   "{contigs.{efasta,fasta,fastb,vecfasta,max.fasta},"
			   "superb,assembly.{fasta,efasta,max.fasta},summary}" ),
		  JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				     "{contigs.vecfasta,superb}" ),
			    FilesIn( sub_dir, "{" + assembly_edits + "}" ) ),

		  "ApplyAssemblyEdits " + pdrs 
		  + ARG( SCAFFOLDS_IN, current_scaffolds)
		  + ARG( EDITS_IN, "\"{" + assembly_edits + "}\"")
		  + ARG( SCAFFOLDS_OUT, current_scaffolds + ".fixed"),
   
		  "Apply edits to the assembly");
    
    make.AddCpRule( sub_dir + "/" + current_scaffolds + ".fixed.rings",
		    sub_dir + "/" + current_scaffolds + ".rings" );
    
    current_scaffolds += ".fixed";
    have_efasta = true;
  }


  if (FIX_LOCAL) {

    // Update alignments if necessary
    if (current_scaffolds != aligned_scaffolds) {
// MakeDepend: dependency AlignReads
      make.AddRule( FilesIn( sub_dir, current_scaffolds + ".readlocs" ),
		    JoinVecs( FilesIn( run_dir, "frag_reads_filt_cpd.{fastb,pairs}",
				       "jump_reads_filt_cpd.{fastb,pairs}"
				       /*
					 ,
					 "long_jump_reads_filt.{fastb,pairs}" 
				       */
				       ),
			      FilesIn( sub_dir, 
				       current_scaffolds + ".{contigs.fastb,superb}" ) ),
		    "AlignReads " + pdrs 
		    + ARG( NUM_THREADS, global_threads ) 
		    + ARG( ASSEMBLY, current_scaffolds ),
		    
		    "Build read locations on the assembly for some of the reads." );

      aligned_scaffolds = current_scaffolds;
    }

// MakeDepend: dependency FixLocal
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".local." +
			   "{contigs.{efasta,fasta,fastb,vecfasta,mapping},"
			   "markup,superb{,.mapping},assembly.{fasta,efasta},summary}" ),
		  JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				     "{contigs.efasta,contigs.fastb,superb,readlocs}" ),
			    FilesIn( run_dir, "ploidy",
				     "frag_reads_filt.{fastb,qualb}",
				     "jump_reads_filt.{fastb,qualb,pairs}",
				     "jump_reads_ec.{fastb,pairs}" ) ),
		  "FixLocal " + pdrs 
		  + ARG( NUM_THREADS, global_threads ) 
            	  + ARG( SCAFFOLDS_IN, current_scaffolds )
                  + ARG( RECYCLE, FL_RECYCLE ),
		  
                  "Locally reassemble and and attempt to correct errors");


    make.AddCpRule( sub_dir + "/" + current_scaffolds + ".local.rings",
		    sub_dir + "/" + current_scaffolds + ".rings" );
    current_scaffolds += ".local";

  }

  if (LONG_READ_POST_PATCH && !BIG_MAP) {
// MakeDepend: dependency PunchTandemHoles
// MakeDepend: dependency LongReadPostPatcher
// MakeDepend: dependency RecoverPunchTandem
   // Note that PunchTandemHoles reads a .markup file so long as FixLocal has been 
   // run.  In a future version, when FixLocal always runs, this should be 
   // explicitly listed as a dependency.

    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".tpunch." +
			   "{edits,tig_ids,{contigs.{efasta,fasta,fastb,vecfasta},"
			   "superb,assembly.{fasta,efasta},summary}}" ),
		  FilesIn( sub_dir, current_scaffolds + "."
		           "{contigs.efasta,superb}" ),
		  "PunchTandemHoles "
		  + ARG(SCAFFOLDS_IN, sub_dir + "/" + current_scaffolds),

		  "Punch holes at tandem repeats preparatory to long read "
                  "patching");

    make.AddCpRule( sub_dir + "/" + current_scaffolds + ".tpunch.rings",
		    sub_dir + "/" + current_scaffolds + ".rings" );
    current_scaffolds += ".tpunch";

    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".longread.edits" ),
		  JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				     "{superb,contigs.{fasta,fastb}}" ),
			    FilesIn( run_dir, "long_reads_orig.fastb" ) ),
		  "LongReadPostPatcher " + pdrs
		  + ARG(NUM_THREADS, global_threads)
		  + ARG(SCAFFOLDS_IN, current_scaffolds)
		  + ARG(SCAFFOLDS_OUT, current_scaffolds + ".longread")
		  + ARG(READS, "long_reads_orig")
                  + ARG(CONTIG_IDS, 
                       "@" + sub_dir + "/" + current_scaffolds + ".tig_ids"),
		  
		  "Patch gaps from PunchTandemHoles using long reads." );

    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".recover."
			   "{superb,contigs.fasta,contigs.fastb,contigs.efasta,"
                           "contigs.vecfasta,assembly.fasta}" ),
		  JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				     "{superb,contigs.fastb,contigs.efasta}" ),
			    FilesIn( sub_dir, current_scaffolds 
                                + ".{edits,longread.edits}" ) ),
		  "RecoverPunchTandem " + 
		  ARG(SCAFFOLDS_IN, sub_dir + "/" + current_scaffolds),
		  
		  "Apply gap patches from PunchTandemHoles/LongReadPostPatcher." );

    make.AddCpRule( sub_dir + "/" + current_scaffolds + ".recover.rings",
		    sub_dir + "/" + current_scaffolds + ".rings" );
    current_scaffolds += ".recover";

  }

  if (KPATCH) {
// MakeDepend: dependency KPatch
    make.AddRule( FilesIn( sub_dir, current_scaffolds + ".kpatch." +
			   "{contigs.{efasta,fasta,fastb,vecfasta,mapping},"
			   "superb{,.mapping},assembly.{fasta,efasta},summary,edits}" ),
		  JoinVecs( FilesIn( sub_dir, current_scaffolds + "."
				     "{contigs.{efasta,fastb},superb}" ),
			    FilesIn( run_dir, "jump_reads_ec.{pairs,fastb}",
                                 ( !BIG_MAP ? "all_reads.unibases.k" + KS
                                 : CLR_uni + ".unibases.k" + ToString(CLR_KOUT)
                                      ) ) ),
		  "KPatch " + pdrs
                  + ARG(K, ( !BIG_MAP ? K : CLR_KOUT ))
                  + ARG(HEAD, ( !BIG_MAP ? "all_reads" : CLR_uni ))
		  + ARG(JUMPS_IN, "jump_reads_ec")
		  + ARG(SCAFFOLDS_IN, current_scaffolds)
		  + ARG(SCAFFOLDS_OUT, current_scaffolds + ".kpatch")
                  + ARG(WRITE_ASSEMBLY, True),

		  "Patch some gaps in the final assembly by walking through "
                  "the unipath graph");

    make.AddCpRule( sub_dir + "/" + current_scaffolds + ".kpatch.rings",
		    sub_dir + "/" + current_scaffolds + ".rings" );
    current_scaffolds += ".kpatch";
  }

  ///////////////////////////////////////////////////////////////////////////
  // Collect together assembly results to a single static name,
  // generating aux assembly files from the superb and efasta.

// MakeDepend: dependency RebuildAssemblyFiles
  make.AddRule( FilesIn( sub_dir, "final.{superb,summary}",
			 "final.contigs.{efasta,fasta,fastb}",
			 "final.assembly.{fasta,efasta}"),
		FilesIn( sub_dir, current_scaffolds + ".{rings,superb,contigs." + 
			 (have_efasta ? "efasta" : "fasta") + "}"),
		"RebuildAssemblyFiles " +
		ARG( IN_HEAD, sub_dir + "/" + current_scaffolds) +
		ARG( OUT_HEAD, sub_dir + "/" + "final"),
		
		"Rebuild final set of assembly files");

  make.AddCpRule( sub_dir + "/final.rings",
		  sub_dir + "/" + current_scaffolds + ".rings" );
  String linear_scaffolds = (have_efasta ? "final.assembly.efasta" : "final.assembly.fasta");


  ///////////////////////////////////////////////////////////////////////////
  // Prepare final assembly for NCBI submission
  // 

// MakeDepend: dependency SubmissionPrep
  make.AddRule( FilesIn( sub_dir, "submission.{agp,"
			 "assembly.{efasta,fasta},superb{,.mapping},"
			 "contigs.{fasta,efasta,mapping}}"),
		FilesIn( sub_dir, "final.{superb,contigs.efasta}"),
		"SubmissionPrep " +
		ARG( HEAD_IN, sub_dir + "/" + "final") +
		ARG( HEAD_OUT, sub_dir + "/" + "submission") +
		ARG( REORDER, SP_REORDER),
		
		"Prepare final assembly for NCBI submission");


  ///////////////////////////////////////////////////////////////////////////
  // Evaluate final assembly
  //

  // Build layout of scaffolds on reference.

  if (evalStandard) {
// MakeDepend: dependency ScaffoldLayout
    make.AddRule( FilesIn( sub_dir, "scaffolds.report" ),
		  JoinVecs( FilesIn( sub_dir, "final.assembly.fasta" ),
			    FilesIn( data_dir, "genome.lookup" ) ),
		  "ScaffoldLayout" + 
		  ARG(NUM_THREADS, global_threads) +
		  ARG(REF_LOOKUP, data_dir + "/genome.lookup") +
		  ARG(ASSEMBLY, sub_dir + "/final.assembly.fasta") +
		  ARG(OUT_HEAD, sub_dir + "/scaffolds") +
		  ARG(TMP_DIR, sub_dir + "/tmp"),

		  "Building scaffold layout on reference" );
  }


  if (evalStandard) {
// MakeDepend: dependency ScaffoldAccuracy
    make.AddOutputRule( FilesIn( sub_dir, "scaffold_accuracy.report" ),
		  JoinVecs( FilesIn( sub_dir, linear_scaffolds ),
			    FilesIn( data_dir, "{genome.fasta,genome.lookup}" ) ),
                        "ScaffoldAccuracy" + 
                        ARG(PERF_STATS, True) +
                        ARG(REFHEAD, data_dir + "/" + SA_REFNAME) +
                        ARG(ASSEMBLY, sub_dir + "/" + linear_scaffolds),
                        
                        "Report on scaffold accuracy" );
  }

  if (evalStandard) {
// MakeDepend: dependency AssemblyCoverage
    make.AddOutputRule( FilesIn( sub_dir, "assembly_coverage.report" ),
                        JoinVecs( FilesIn( sub_dir, linear_scaffolds ),
                                  FilesIn( data_dir, "genome.{lookup,fasta}" ) ),
                        "AssemblyCoverage" + 
                        ARG(NUM_THREADS, global_threads) +
                        ARG(REF, data_dir + "/genome.fasta") +
                        ARG(ASSEMBLY, sub_dir + "/" + linear_scaffolds) +
                        ARG(HEAD, sub_dir + "/AssemblyCoverage"),
                        
                        "Report on assembly coverage" );
  }

  if (evalStandard) {
// MakeDepend: dependency AssemblyAccuracy
    make.AddOutputRule( FilesIn( sub_dir, "assembly_accuracy.report" ),
                        JoinVecs( FilesIn( sub_dir, linear_scaffolds ),
                                  FilesIn( data_dir, AA_REFNAME + ".fasta" ) ),
                        "AssemblyAccuracy" +
                        ARG(NUM_THREADS, global_threads) +
                        ARG(REF, data_dir + "/" + AA_REFNAME + ".fasta") +
                        ARG(CHUNK_SIZE, AA_CHUNK_SIZE) +
                        ARG(ASSEMBLY, sub_dir + "/" + linear_scaffolds),
                        
                        "Report on assembly accuracy" );
  }

  if (evalStandard) {
// MakeDepend: dependency EvaluateGaps
    make.AddOutputRule( FilesIn( sub_dir, "assembly_gaps.report" ),
                        JoinVecs( FilesIn( sub_dir, 
                             "final.{contigs.fastb,superb}" ),
                                  FilesIn( data_dir, "genome.lookup" ) ),
                        "EvaluateGaps " + pdrs + 
			ARG(NUM_THREADS, global_threads) +
			ARG(SCAFFOLDS_IN,"final"),
                        
                        "Report on assembly gaps" );
  }

  // temporarily add gap evaluation on linear_scaffolds0
  if (evalStandard && !BIG_MAP) {
// MakeDepend: dependency EvaluateGaps
    make.AddRule( FilesIn( sub_dir, "linear_scaffolds0.contigs.fastb" ),
		  FilesIn( sub_dir, "linear_scaffolds0.contigs.fasta" ),
		  "Fasta2Fastb "  
		  " NAMES=False" +
		  ARG(IN, sub_dir + "/linear_scaffolds0.contigs.fasta"),

		  "Convert linear_scaffolds0.contigs.fasta to fastb" );

    make.AddOutputRule( FilesIn( sub_dir, "initial_assembly_gaps.report" ),
                        JoinVecs( FilesIn( sub_dir, 
                             "linear_scaffolds0.{contigs.fastb,superb}" ),
                                  FilesIn( data_dir, "genome.lookup" ) ),
                        "EvaluateGaps " + pdrs + 
			ARG(NUM_THREADS, global_threads) +
			ARG(SCAFFOLDS_IN,"linear_scaffolds0"),
                        
                        "Report on assembly gaps" );
  }
  if (evalStandard) {
// MakeDepend: dependency EvalScaffolds
    make.AddRule( FilesIn( sub_dir, "EvalScaffolds/main.log" ),
		  JoinVecs( FilesIn( sub_dir, "final.{superb,contigs.fastb}" ),
			    FilesIn( data_dir, "genome.{lookup,fastb}" ) ),
		  "EvalScaffolds" + 
		  ARG(NUM_THREADS, global_threads) +
		  ARG(LOOKUP, data_dir + "/genome.lookup") +
		  ARG(SCAFFOLDS, sub_dir + "/final") + 
		  ARG(OUT_DIR, sub_dir + "/EvalScaffolds"),

		  "Evaluate scaffolds using reference" );
  }

  // General reports

// MakeDepend: dependency AllPathsReport
  make.AddOutputRule( FilesIn( sub_dir, "assembly_stats.report" ),
		      FilesIn( sub_dir, "final.{superb,contigs.efasta}" ),
		      "AllPathsReport " + pdrs +
		      ARG(ASSEMBLY, "final"),
			
		      "Generate basic stats - without a reference" );


// MakeDepend: dependency LibCoverage
  make.AddOutputRule( FilesIn( sub_dir, "library_coverage.report" ),
		      JoinVecs( FilesIn( sub_dir, aligned_scaffolds + 
					 ".{superb,readlocs,contigs.fastb}" ),
				FilesIn( run_dir, "" ) ),
		      "LibCoverage " + pdrs +
		      ARG(ASSEMBLY, aligned_scaffolds) +
		      ARG(NUM_THREADS, global_threads),
		      
		      "Generate basic library coverage stats - without a reference" );


  //
  //  End of Rule Definitions
  //
  /////////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////
  // Specify targets to build
  //

  vec< filename_t > targetList;

  // Add specified target files
  targetList.append( FilesIn( run_dir, TARGETS_RUN ));
  targetList.append( FilesIn( sub_dir, TARGETS_SUBDIR ) );
  targetList.append( FilesIn( data_dir, TARGETS_DATA ) );
  targetList.append( FilesIn( run_dir, TARGETS_REF ));

  // Add pseudo targets - duplication here doesn't matter

  // Targets for NCBI submission
  if ( Member(TARGETS, String("submission"))  ) 
    TARGETS.push_back("standard");


  // Targets for full eval case (full evaluation + standard)
  if ( Member(TARGETS, String("full_eval"))  ) 
    TARGETS.push_back("standard");


  // Special targets 
  
  // filled_fragments - Generate filled fragments only
  if ( Member(TARGETS, String("filled_fragments"))  ) 
    targetList.append(FilesIn(run_dir, "filled_reads.{fastb,pairs}") );

  // frag_unipaths - Generate fragment unipaths only
  if ( Member(TARGETS, String("frag_unipaths"))  ) 
    targetList.append(FilesIn(run_dir, "extended40.shaved.unipaths.k96") );

  // // LocalizeReads - Run up to, but not including,  LocalizeReadsLG
  // if ( Member(TARGETS, String("LocalizeReads")) ) 
  //   targetList.append(LocalizeReadsDeps);

  // submission - Prepare for NCBI submission
  if ( Member(TARGETS, String("submission"))  ) 
    targetList.append(FilesIn(sub_dir, "submission.assembly.efasta") );

  // General targets
    
  // Standard default targets 
  if ( Member(TARGETS, String("standard"))  ) {

    // Actual assembly files
    targetList.append(FilesIn(sub_dir, "{"
			      "final.{assembly.fasta,superb,rings},"
			      "final.contigs.{fasta,fastb}"
			      "}" ) );
    if (have_efasta) 
      targetList.append(FilesIn(sub_dir, "final.{assembly,contigs}.efasta" ) );
    
    // basic evaluation
    if (evalBasic) {
      targetList.append(FilesIn(sub_dir, "{"
				"assembly_stats.report,"
				"library_coverage.report,"
				"}" ) );
      //      if (BIG_MAP)
	targetList.append(FilesIn(run_dir, "jump_reads_ec.distribs")); 
    }

    // standard evaluation
    if ( evalStandard ) {
      targetList.append(FilesIn(sub_dir,"{" 
				"scaffold_accuracy.report,"
                                "assembly_accuracy.report,"
                                "assembly_gaps.report,"
                                "assembly_coverage.report,"
                                "scaffolds.report,"
				"EvalScaffolds/main.log,"
                                + ToString( !BIG_MAP ? "seeds.map" : "" ) + "," 
				"}" ) );
      if (!BIG_MAP) 
     	targetList.append(FilesIn(sub_dir, "initial_assembly_gaps.report")); 
    }
  }

  // Targets for full evaluation
  if ( Member(TARGETS, String("full_eval"))  ) {
    if ( evalStandard ) {
      targetList.append(FilesIn(run_dir,"{all_reads.unipaths.kN.{UnipathEval.out,placements}}" ));
      if ( !BIG_MAP )
           targetList.append(FilesIn(sub_dir,"EvalUnipathLinkGraphs.out") );
      //      targetList.append( FilesIn( sub_dir, aligned_scaffolds + ".u2c" ) );
    }
  }


  // Additional targets
  if ( VALIDATE_INPUTS )
    targetList.append(FilesIn( data_dir, "ValidateAllPathsInputs.report" ));


  if (UNIPATH_BABY_STATS)
    targetList.append(FilesIn(run_dir, "UnipathBabyStats.out" ));


  if ( FORCE_TARGETS && !Member(TARGETS , String("all")))
    make.AddForceTargets( targetList );
  
  make.AddForceTargetsOf( FORCE_TARGETS_OF );
  
  if ( RESTORE_TARGETS_OF.nonempty() )
    make.RestoreTargets();

  if (Member(TARGETS, String("all")))
    // special pseudo target "all" - makes all know targets
    targetList = MkVec(ToString("all"));
  
  // Set targets to make
  make.SetTargets( targetList, vec< filename_t >() );
  
  // Write dependancy graph - seems to be broken right now.
  // make.WriteDependencyGraph( make_mgr_root + ".dot" );

  if ( VIEW_PIPELINE_AND_QUIT ) {
    make.WritePipelineDescrAsHTML( log_dir + "/pipeline" );
    
    cout << "Pipeline written out, quitting..." << endl;
    exit( 0 );
  }

  cout << endl << Date() << ": Starting ALLPATHS-LG Pipeline." << endl << endl; 

  // Test rules (or rule) and stop
  if (TEST_RULE != "") 
    return make.TestRules(ref_dir, TEST_RULE, TEST_MODE, TEST_ERASE );
  
  // Run the pipeline 
  int status = make.RunMake( MAXPAR );

  if (status == 1)  // nothing to be done
    cout << endl << "Requested targets are up to date - nothing to be done." << endl;

  else if (DRY_RUN == False) {  // Generate reports

// MakeDepend: dependency install_scripts
// MakeDepend: dependency MemMonitor
    String parseMMCmd = "ParseMemMonitorOutput.pl";

    if (MEMORY_DIAGNOSTICS && IsCommandInPath(parseMMCmd)) {
      cout << endl << Date() << ": Computing runtime statistics." << endl << endl; 
      SystemSucceed( parseMMCmd + " SUBDIR=" + sub_dir);
    }

// MakeDepend: dependency install_scripts
    String CARCmd ="CompileAssemblyReport.pl";
    
    if (IsCommandInPath(CARCmd)) {
      cout << endl << Date() << ": Compiling assembly report." << endl; 
      SystemSucceed( CARCmd + " NH=True SUBDIR=" + sub_dir);
    }
  }
  
  cout << endl << Date() << " : ALLPATHS-LG Pipeline Finished." << endl << endl;

  cout << "Run directory: " << run_dir << endl;
  cout << "Log directory: " << log_dir << endl;

  if ( status == 2 ) {
    cout << endl << " *** Make encountered an error, see above for error messages. ***"
	 << endl << endl;;
    // Return error code if pipeline failed during make stage
    return 1;
  }  
}
