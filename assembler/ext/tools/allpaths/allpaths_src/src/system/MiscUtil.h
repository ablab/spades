///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Header file: MiscUtil.h

   Miscellaneous useful utilities -- put into this file to avoid creating
   unnecessary dependencies.
   
   @file
*/

#ifndef __INCLUDE_system_MiscUtil_h
#define __INCLUDE_system_MiscUtil_h

#include <set>
#include <map>
#include "Vec.h"
#include "system/Types.h"
#include "system/System.h"
#include "system/HTMLUtils.h"
#include "CommonSemanticTypes.h"
#include "graph/Digraph.h"

// Semantic type: shellcmd_t
// The name of a shell command.
SemanticType( String, shellcmd_t );

/**
    FuncDecl: NeedToMake

    Given a set of target files and a set of source files,
    returns True if any of the target files is missing or
    or older than one of the source files.

    Halts with a fatal error if any of the source files
    are missing.
*/
Bool NeedToMake( const vec< filename_t >& targetFiles,
		 const vec< filename_t >& sourceFiles );


/**
   FuncDecl: Make

   Make target files from source files using the specified
   command.  The command is run only if needed,
   unless runOnlyIfNeeded==False.
*/
void Make( const vec< filename_t >& targetFiles,
	   const vec< filename_t >& sourceFiles,
	   shellcmd_t command,
	   Bool runOnlyIfNeeded = True );

/**
   FuncDecl: Make

   Make target files from source files using the specified
   command.  The command is run only if needed,
   unless runOnlyIfNeeded==False.
*/
void Make( const vec< filename_t >& targetFiles,
	   const vec< filename_t >& sourceFiles,
	   const vec< shellcmd_t >& commands,
	   Bool runOnlyIfNeeded = True );

void SetMakeOnlyIfNeeded( Bool );


/**
   FuncDecl: PrependDir

   Prepend a String to a vec of Strings -- for example, prepends a dirname to each filename.
*/
vec< filename_t > PrependDir( const filename_t& dir, const vec< filename_t >& files );


/**
   FuncDecl: AbsToRelPaths

   Convert two absolute paths to a relative path.
*/
void AbsToRelPaths(const String& abs_source, const String& abs_dest, String& rel_source);

/**
   FuncDecl: FilesIn

   Make a vector of filenames in the given dir.
 */
vec< filename_t > FilesIn( const dirname_t& dir, const filename_t& fn1 );
vec< filename_t > FilesIn( const dirname_t& dir, const filename_t& fn1, const filename_t& fn2 );
vec< filename_t > FilesIn( const dirname_t& dir, const filename_t& fn1, const filename_t& fn2,
			   const filename_t& fn3);
vec< filename_t > FilesIn( const dirname_t& dir, const filename_t& fn1, const filename_t& fn2,
			   const filename_t& fn3, const filename_t& fn4 );
vec< filename_t > FilesIn( const dirname_t& dir, const filename_t& fn1, const filename_t& fn2,
			   const filename_t& fn3, const filename_t& fn4, const filename_t& fn5 );
vec< filename_t > FilesIn( const dirname_t& dir, const filename_t& fn1, const filename_t& fn2,
			   const filename_t& fn3, const filename_t& fn4, const filename_t& fn5,
			   const filename_t& fn6 );
vec< filename_t > FilesIn( const dirname_t& dir, const filename_t& fn1, const filename_t& fn2,
			   const filename_t& fn3, const filename_t& fn4, const filename_t& fn5,
			   const filename_t& fn6, const filename_t& fn7 );
vec< filename_t > FilesIn( const dirname_t& dir, const filename_t& fn1, const filename_t& fn2,
			   const filename_t& fn3, const filename_t& fn4, const filename_t& fn5,
			   const filename_t& fn6, const filename_t& fn7, const filename_t& fn8 );
vec< filename_t > FilesIn( const dirname_t& dir, const filename_t& fn1, const filename_t& fn2,
			   const filename_t& fn3, const filename_t& fn4, const filename_t& fn5,
			   const filename_t& fn6, const filename_t& fn7, const filename_t& fn8,
			   const filename_t& fn9 );

void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1 );
void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2 );
void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2, filename_t fn3 );
void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2, filename_t fn3, filename_t fn4 );
void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2, filename_t fn3, filename_t fn4, filename_t fn5 );
void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2, filename_t fn3, filename_t fn4, filename_t fn5, filename_t fn6 );
void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2, filename_t fn3, filename_t fn4, filename_t fn5, filename_t fn6,
		filename_t fn7 );
void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2, filename_t fn3, filename_t fn4, filename_t fn5, filename_t fn6,
		filename_t fn7, filename_t fn8 );
void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2, filename_t fn3, filename_t fn4, filename_t fn5, filename_t fn6,
		filename_t fn7, filename_t fn8, filename_t fn9 );


/**
    Class: MakeRule

    Represents one rule in a MakeMgr.  A rule specifies how to make a set of target data files
    from a set of source data files, by running a command.  With each command there is also
    an associated comment describing what the command does.

    This class is internal to MakeMgr.  As a user, you never create instances of MakeRule;
    instead you call methods of MakeMgr -- MakeMgr::AddRule(), MakeMgr::AddCpRule() etc --
    to add rules.
 */
class MakeRule {
 public:
  MakeRule() { }
  MakeRule( const vec< filename_t >& targetFiles_, const vec< filename_t >& sourceFiles_,
	    shellcmd_t command_, filename_t saveOutputTo_, String comment_,
	    String shortRuleName_, String mediumRuleName_ );

  const vec< filename_t >& Targets() const { return targetFiles; }
  const vec< filename_t >& Sources() const { return sourceFiles; }

  /// Method: Command
  /// Returns the shell command that needs to be run to make the targets from the sources.
  const shellcmd_t& Command() const { return command; }

  /// Method: SaveOutputTo
  /// Returns the name of the (optional) file to which the output of this rule's command
  /// is to be logged.
  filename_t SaveOutputTo() const { return saveOutputTo; }

  shellcmd_t CommandName() const { return WhiteSpaceFree( Command().SafeBefore( " " ) ); }

  String Comment() const { return comment; }
  String CommentOrCommand() const { return comment.nonempty() ? comment : command; }
  String ShortRuleName() const { return shortRuleName; }
  String MediumRuleName() const { return mediumRuleName; }

  /// Method: PrimaryTarget
  /// Returns the primary target of the command; we arbitrarily choose one of the command's targets
  /// -- here the first target -- as the primary target.  See the documentation of the
  /// targetFiles field for explanation.
  filename_t PrimaryTarget() const { return targetFiles.front(); }
  
 private:
  /// Field: targetFiles
  /// The files written by this command.
  /// After the command is run, we check that it has in fact created
  /// all these files.
  /// 
  /// One target of the rule is chosen as the _primary target_.  It is returned
  /// by the PrimaryTarget() method.  Dependencies on any targets of this rule
  /// are changed in the written-out makefile into dependencies on the primary target.
  /// This lets us correctly handle the case where a command writes out several targets,
  /// and parallel make is used.   If we wrote out a separate rule for each target,
  /// make would start several copies of the command, one for making each target,
  /// as required by downstream targets.   This would both waste resources and
  /// lead to errors as several instances of the command compete to write out
  /// the same set of targets.   Designating one target as the primary target,
  /// and changing all dependencies on targets of this command into dependencies
  /// on the primary target, solves this problem.
  vec< filename_t > targetFiles;
  
  /// Field: sourceFiles
  /// The files read by this command.
  vec< filename_t > sourceFiles;
  
  /// Field: command
  /// The shell command to invoke, to create the targets from the sources.
  /// Usually the command is either an Arachne/ALLPATHS program or the cp command.
  shellcmd_t command;

  /// Field: saveOutputTo
  /// Optionally, if we are saving/logging a copy of all output of the command to a file,
  /// the name of that file.
  filename_t saveOutputTo;

  /// Field: comment
  /// A user-readable description of what the command does.  It is printed before the command
  /// is run, and also shown in the pipeline visualization in MakeMgr::WritePipelineDescrAsHTML().
  String comment;

  /// Field: shortRuleName
  /// A short string that uniquely identifies this rule within a MakeMgr.   Normally composed of
  /// the first letters of words used in the command (e.g. LocalizeReads becomes LR).  Used e.g. for prefixing
  /// printed lines relating to this command when the makefile is run (this helps when several
  /// commands are run in parallel and all write to the same console -- this helps tell which line
  /// relates to which command).
  String shortRuleName;

  /// Field: mediumRuleName
  /// A medium-length string that uniquely identifies this rule within a MakeMgr.
  /// Normally composed of the command name, with a number added at the end for disambiguation
  /// if needed (if there are several rules with the same command in a MakeMgr).
  /// Used in various output settings, for example for drawing the pipeline as a graph
  /// in MakeMgr::WritePipelineDescrAsHTML().
  String mediumRuleName;
};  // class MakeRule

/**
   Class: MakeMgr

   A Makefile facility specialized not for compiling a codebase but for running a research analysis pipeline.
   You programmatically (by calling methods of MakeMgr) specify a set of datafiles and set of programs,
   where each program reads some datafiles (sources) and produces others (targets).
   MakeMgr then creates a makefile representing all the dependencies, and
   runs make.

   Thus, you get makefile-like functionality (only targets that need to be rebuilt are in fact rebuilt,
   and independent targets can be made in parallel), but with various additional benefits.

   The program RunAllPaths is an example of MakeMgr usage.

   Functionality provided by MakeMgr:

     - it catches errors early: makes sure the entire pipeline is consistent
     (that each target is either declared pre-existing or we're told how to make it),
     that there are no circular dependencies, etc.   These checks are made before
     anything is run.  This matters because some runs can take days, re-running is expensive,
     so it's important to know early if the run is guaranteed to run into problems.

     Additional checks are easy to add.  For example, a useful check would be to prevent
     conflicts between several people trying to create the same files.  When a pipeline run
     is started we know all the files that will eventually be written, so we can lock them
     all until the pipeline finishes, and prevent others from erroneously starting pipeline
     runs that will write the same files.

     - it adds run-time error checking -- for example after each command is run, we check
     that is has in fact created all the files it was supposed to create, and abort immediately
     if it didn't.  This catches errors long before these missing files would cause crashes
     later in the pipeline.

     - if a program crashes or is interrupted, and leaves partially
     written files, these files will be recognized as incomplete
     and automatically deleted.

     - the pipeline can be graphically visualized as an HTML page.
     This helps in planning changes/extensions to the pipeline.

     - for each target file, we automatically record various useful metadata
     about how the file was made.  we record the command that was run,
     the sources and targets of the command,
     and the output of that command.  for a file D/F (file named F in directory D),
     the metadata files are named D/makeinfo/F.*  In the future, metadata might
     include how long the command took to run, whether the file will be overwritten
     by a currently running pipeline run, md5 checksums of source files on which
     the file depends, info about the file's format, etc.
     All this helps reconstruct how the file was made, even if the pipeline
     changes in the future.

     - MakeMgr supports a dynamically changing research pipeline.  as noted above,
     for each target file we record the command used to make that file, and the
     source file dependencies.
     if the command line changes, or new source file dependencies are added,
     the file is forced to be remade.  at the same time, it is recognized that
     some command-line parameter changes do not affect the output files (e.g.
     debugging options, or options specifying how much parallelization is to be used
     during the computation).   so MakeMgr has a facility to ignore changes
     to certain parameters (see AddIgnoreChangesToCommandArgs() ).

     In an experimental research pipeline, there are often two or more alternate ways
     of computing a given file.  The metadata keeps track of how the file was last made
     (i.e. by running what command line), which lets us correctly keep files up-to-date in the
     face of changes to the pipeline.

   Compared to writing your own makefiles, MakeMgr adds the following functionality:

     - GNU make does not recognize that a command may write multiple targets.
     You can create separate rules for each target of a command, but this will break
     if you use the automatic parallelization function of make (-j6).  Make may start
     several copies of the command to make each of its targets, rather than recognizing
     that only one copy of the command needs to be invoke to make all the targets.
     The several copies of the command will waste resources and, worse, will all be writing
     the same set of files, which is likely to lead to problems.

     By contrast, MakeMgr explicitly knows that your command makes several targets at once.
     It adjusts the makefile rules to make sure this is handled correctly even in the face
     of parallel make.  Specifically, it picks one target of the command as the "canonical"
     target; any dependencies by other files on any target of the command become dependencies
     on the canonical target.  Thus, the canonical target becomes a proxy for whether the entire
     set of the command's targets is up-to-date.

     - There is much boilerplate involved in writing a makefile.  Specifically, in the context
     of Arachne/ALLPATHS, there are often many files with similar names.  MakeMgr lets you specify
     sets of similarly-named files using the ParseStringSet() syntax.

     - The makefile language is much more general than what is needed for the research pipeline,
     and is tailored for describing how to compile a codebase rather than in encoding the knowledge
     of how to compute various derived data files from primary data files.
     By providing only the functionality needed for a research pipeline, we get correspondingly better error-checking.
     On the other hand, using C++ to construct the rules programmatically gives us the flexibility
     in constructing command lines and in allowing for alternate ways of making a given file.

     It also lets us recognize and handle some common special cases -- for example, special handling
     for the simple class of rules that obtain a target file by simply copying a source file.

   Usage notes:

   See program RunAllPaths for examples of usage.

   In general, there are three stages of using a MakeMgr.
   In the first stage, you call various configuration methods, e.g. to define
   abbreviations for commonly used subdirectories.
   Then you call PrepareForRules().
   In the second stage, you add rules by calling AddRule() or AddCpRule().
   Each rule specifies how to make a set of target files from a set of source files,
   using a specified command.  Both the source and target file sets, and the
   command, are constructed programmatically, which gives useful flexibility.
   In the third stage, you call RunMake() to construct the makefile based on the rules,
   and to invoke make to run the commands.

   Implementation notes:

   Vectors of filenames are small, so sometimes we pass them around inefficiently --
   using copy-by-value to copy the entire vector -- for the sake of more compact code.
   Also, internally we always pass around and represent filenames as absolute full
   paths, even though they are insanely long if you print them.
   For printing, you can pass the filenames through ApplySubsts() which will substitute
   the common directory names with abbreviations specified by the user (e.g. $(RUN)/filename).

   To-do:

      - automatically add dependencies of targets on the C++ program binaries
      used to create the targets
      - and run make to rebuild the C++ program binaries on which we depend.
      - automatically add dependencies on the Arachne program parameters
      (independent of program order); have an option to just warn if
      program params have changed.
*/
class MakeMgr {
 public:

  /// MakeMgr constructor
  /// The makefile will be created and run in /tmp.  To save a copy of this
  /// makefile elsewhere, pass an appropriate filename to the constructor.
  MakeMgr(filename_t out = "mytest.make") {
    outFilename = out;
    preparedForRules = False;
    lnForCp = False;
    ignoreMetaData = False;
    dryRun = False;
    continueAfterError = False;
    diagnostics = False;
    ignoreInitialPartialOutputs = False;
    logDirectory = "";
    makeCommand = "make";
  }

  // Section: Configuring the MakeMgr
  //
  // Methods for configuring the MakeMgr before you start
  // adding rules.

  /// Method: AddDirAbbrev
  /// Specify an abbreviated name for a directory.
  /// For example, you might abbreviate a run dir
  /// /wga/dev/WGAdata/projects/ALLPATHS/E.coli/jaffe16
  /// as RUN.  Then references to file F in that directory
  /// will be written out as $(RUN)/F whenever the filename
  /// is displayed for human consumption.
  /// Sample usage:
  /// >   make.AddDirAbbrev( ref_dir, "REF" );
  /// >   make.AddDirAbbrev( data_dir, "DATA" );
  /// Note that if DATA is a subdir of REF, files in the data dir
  /// will be correctly abbreviated as $(DATA)/F and not as
  /// $(REF)/something/F.
  /// Note also that MakeMgr knows nothing about our particular
  /// directory organization used in ALLPATHS (run dir, data dir, ...).
  /// We simply tell is the set of common directories that we use,
  /// and their abbreviations.
  void AddDirAbbrev( dirname_t dir, String abbr ) {
    ForceAssert( !preparedForRules );
    dirAbbrevs[ dir ] = abbr;
  }

  /// Method: AddSubst
  /// Add a replacement to be applied to all filenames.
  /// Occurrences of the 'from' string will be replaced by the 'to'
  /// string.
  /// Sample usage:
  /// >  make.AddSubst( ".kN", ".k" + ToString( K ) );
  /// where K is the kmer size used. 
  /// This allows more compact programmatic construction of filenames --
  /// you just write "myfile.kN.ext" rather than "mfile.k" + ToString( K ) + ".ext".
  void AddSubst( String from, String to ) {
    ForceAssert( !preparedForRules );
    substs[ from ] = to;
  }

  /**
     Method: AddPreExistingFile

     Add the name of a file which we can expect to exist when we start the make,
     and so which we do not need to know how to build.  Dependencies on this
     file are therefore not errors even if we don't have a rule for making
     this file.

     If we do have a rule for making this file, the rule will be ignored and the
     file will still be treated as pre-existing.
  */
  void AddPreExistingFile( filename_t fn ) {
    ForceAssert( !preparedForRules );
    preExistingFiles.insert( fn );
  }

  /// Method: AddPreExistingFiles
  /// Add multiple preexisting files; see AddPreExistingFile().
  void AddPreExistingFiles( vec< filename_t > fnames ) {
    ForceAssert( !preparedForRules );
    fnames = ExpandFileNames( ApplySubsts( fnames ) );
    preExistingFiles.insert( fnames.begin(), fnames.end() );
  }

  /**
      Method: AddDontUpdateTargetsOf

      Specify that the targets of the given command(s) should not be updated,
      even if they're out-of-date.  Useful when remaking the targets involves
      running time-consuming commands, and you just want to run the pipeline from
      a given point forward and experiment with files that are downstream
      from the targets of the specified commands.

      The rules are specified by giving the medium command name; if the command
      is named CMD and there is just one such command (i.e. just one invocation
      of this program) in the MakeMgr then the medium rule name is CMD,
      otherwise it is CMD-1 or CMD-2 etc -- which one it is you can determine
      by looking at the HTML representation of the pipeline (written out by
      WritePipelineDescrAsHTML(); see the RunAllPaths command-line option
      VIEW_PIPELINE_AND_QUIT).
   */
  void AddDontUpdateTargetsOf( const vec< shellcmd_t >& commandNames ) {
    ForceAssert( !preparedForRules );
    dontUpdateTargetsOf.insert( commandNames.begin(), commandNames.end() );
  }

  void AddRestoreTargetsOf( const vec< shellcmd_t >& commandNames ) {
    ForceAssert( !preparedForRules );
    restoreTargetsOf.insert( commandNames.begin(), commandNames.end() );
  }

  /**
      Method: AddForceTargetsOf

      Specify commands whose targets will be force to be remade, even if they are
      up-to-date.  See the RunAllPaths command-line option FORCE_TARGETS_OF.
      Commonly you'd use this if you changed the program run by this rule,
      in some functionally important way.

      We could of course make the rule's target files dependent on the
      date of the rule's executable, and force them to be remade if
      the executable changes.  But recompilations of executables is
      very conservative, and often an executable gets recompiled even
      though functionally it has not changed.   So forcing us to rerun
      a command if its executable has changed would lead to many unneeded
      reruns of (possibly very long-running) commands.  So we don't do that,
      but the flipside is that when a program _is_ changed in a functionally
      significant way, you have to manually force its targets to be remade.
   */
  void AddForceTargetsOf( const vec< shellcmd_t >& commandNames ) {
    forceTargetsOf.insert( commandNames.begin(), commandNames.end() );
  }

  /**
      Method: AddForceTargets

      Force the specified target files to be remade.  See comments for
      the method AddForceTargetsOf().
   */
  void AddForceTargets( const vec< filename_t >& targetsToForce ) {
    forceTheseTargets.append( targetsToForce );
  }

  /**
     Method: IgnoreMissingMetadataFor

     For the files listed, if the metadata ( makeinfo/FILENAME.cmd ) about how
     the file was made is missing, do not force a rebuild of this file.
     The metadata may be missing if you manually directly ran the command,
     rather than invoking it through RunAllPaths.   Or if the command used to
     write out the target files but they were not explicitly listed as targets
     of the command in the corresponding MakeMgr rule.
     In either case, this is a hack that bypasses the normal update mechanisms --
     used it to save unnecessary reruns of commands, but use at your own risk.
  */
  void IgnoreMissingMetadataFor( vec< filename_t > ignoreMissingMetadataFor_ ) {
    ignoreMissingMetadataFor.append( ExpandFileNames( ApplySubsts( ignoreMissingMetadataFor_ ) ) );
  }


  void SetLnForCp( Bool lnForCp_ ) {
    ForceAssert( !preparedForRules );
    lnForCp = lnForCp_;
  }

  void SetIgnoreMetaData( Bool ignoreMetaData_ ) {
    ForceAssert( !preparedForRules );
    ignoreMetaData = ignoreMetaData_;
  }

  void SetDryRun( Bool dryRun_ ) {
    ForceAssert( !preparedForRules );
    dryRun = dryRun_;
  }

  void SetMakeCommand( const String& makeCommand_ ) {
    ForceAssert( !preparedForRules );
    makeCommand = makeCommand_;
  }

  void SetContinueAfterError( Bool continueAfterError_ ) {
    ForceAssert( !preparedForRules );
    continueAfterError = continueAfterError_;
  }

  void SetDiagnostics( Bool diagnostics_ ) {
    ForceAssert( !preparedForRules );
    diagnostics = diagnostics_;
  }

  void SetIgnoreInitialPartialOutputs( Bool ignoreInitialPartialOutputs_ ) {
    ForceAssert( !preparedForRules );
    ignoreInitialPartialOutputs = ignoreInitialPartialOutputs_;
  }

  void SetLogDirectory( const String& logDirectory_ ) {
    ForceAssert( !preparedForRules );
    logDirectory = logDirectory_;
  }

  void AddIgnoreChangesToCommandArgs( String commandArgs );

  /// Method: AddSpecialTarget
  /// Add a target that represents a specified set of files to make.
  /// You might use it to create predefined targets for making some common
  /// groups of files, as is commonly done in makefiles.
  /// The target is "special" because it is not the name of any real file,
  /// but a logical name for a set of target files.
  void AddSpecialTarget( String specialTarget, vec< filename_t > filesInTarget );

  ////////////////////////////////////

  // Section: Adding rule to MakeMgr
  //
  // The following methods are invoked after you're done configuring the MakeMgr.

  /**
      Method: PrepareForRules

      Tell MakeMgr that we have finished calling all the configuration routines on it,
      and will now be adding the actual rules.

      After this, configuration routines may no longer be called.
   */
  void PrepareForRules();

  /**
     MethodDecl: MakeRule
  
     Define a rule for how to make a set of target files
     from a set of source files by running a given command.
  */
  void AddRule( vec< filename_t > targetFiles,
		vec< filename_t > sourceFiles,
		shellcmd_t command,
		String comment = "",
		filename_t saveOutputTo = ""
		);

  void AddRule( vec< filename_t > targetFiles,
		vec< filename_t > sourceFiles,
		shellcmd_t command,
		String comment,
		vec< filename_t > saveOutputTo ) {
    ForceAssert( saveOutputTo.solo() );
    AddRule( targetFiles, sourceFiles, command, comment, ApplySubsts( saveOutputTo.front() ) );
  }

  void AddOutputRule( vec< filename_t > targetFiles,
		      vec< filename_t > sourceFiles,
		      shellcmd_t command,
		      String comment = "" ) {
    AddRule( targetFiles, sourceFiles, command, comment,
	     ApplySubsts( targetFiles.front() ) );
  }
  

  void AddCpRule( filename_t target, filename_t source );
  void AddCpRule( dirname_t targetDir, dirname_t sourceDir, const vec< filename_t >& fnames );
  void AddCpRule( dirname_t targetDir, filename_t targetStart, 
		  dirname_t sourceDir, filename_t sourceStart,
		  const filename_t  commonEnd );
  void AddCpRule( dirname_t targetDir, dirname_t sourceDir, filename_t fn1 );

  void SetTargets( vec< filename_t > defaultTargets,
		   vec< filename_t > excludeTargets );

  // MethodDecl: RunMake
  // Run the make process. Returns modified make exit code:
  // 0 = make sucessful; 1 = targets up to date; 2 = error occured
  int RunMake( int MAX_PARALLEL = 4 );
  

  int TestRules( String base_name, String ref_name, String mode, bool erase);
  
  void WriteDependencyGraph( filename_t dotFileName ) const {
    WriteDependencyGraph( dotFileName, True, NULL, NULL, NULL, NULL, "Deps" );
  }

  void RestoreTargets();

  /**
     MethodDecl: WritePipelineDescrAsHTML
     
     Create a webpage that describes the entire pipeline.
  */
  void WritePipelineDescrAsHTML( dirname_t dir ) const;
  
  
  filename_t ApplySubsts( filename_t fn ) const;
  vec< filename_t > ApplySubsts( const vec< filename_t >& fnames ) const;

  
  vec< filename_t > ExpandFileNames( const vec< filename_t >& fnames ) const;

  //////////////////////////////////////////////////////////////////////
  
 private:

  /// Field: rules
  /// The rules; see MakeRule.  Each rule specifies how to make a set of target files
  /// from a set of source files by invoking a given command.
  vec< MakeRule > rules;

  // Type: dir2abbrev_t
  // A map of abbreviated names for commonly-used directories.
  typedef map< dirname_t, String > dir2abbrev_t;

  /// Field: dirAbbrevs
  /// A map of abbreviated names for commonly-used directories.
  /// Most data files reside in a handful of directories,
  /// so for display purposes it is useful to abbreviate long
  /// directory names using these abbreviations.
  dir2abbrev_t dirAbbrevs;

  /// Field: abbreviatedDirsRevByLen
  /// The names of directories for which the user has defined abbreviated names,
  /// sorted by length.  Often several abbreviated directories are subdirs of
  /// one another, so this helps us find the right abbreviation for a given filename.
  vec< dirname_t > abbreviatedDirsRevByLen;

  /// Field: substs
  /// Substitutions to be applied to each file name, e.g.
  /// .kN to .k20 .
  map< String, String > substs;

  // Type: rule_id_t
  // Identifier of a rule within a MakeMgr; the index of that MakeRule in the 'rules' field
  // of the MakeMgr.
  typedef int rule_id_t;

  /**
     Field: target2rule

     For each target, the id of the rule showing how to make
     that target.
  */
  map< filename_t, rule_id_t > target2rule;

  /// Logging directory
  /// directory where the logs generated by each pipeline run will be written
  String logDirectory;

  /// Field: outFilename
  /// The filename to which the generated makefile should be copied.
  filename_t outFilename;

  /// Field: preparedForRules
  /// Whether PrepareForRules() has been called, indicating that all preliminary
  /// configuration of the MakeMgr is complete and that we are now adding rules to it.
  Bool preparedForRules;

  /// Field: preExistingFiles
  /// These files are assumed to already exist.  We check that they do in fact exist
  /// before we start running the pipeline.
  set< filename_t > preExistingFiles;

  /// Field: shortRuleNames
  /// For constructing unambiguous short rule names, unique within this MakeMgr: the short rule
  /// names already used.
  set< String > shortRuleNames;
  
  /// Field: mediumRuleNames
  /// For constructing unambiguous medium rule names, unique within this MakeMgr: the medium rule
  /// names already used.
  set< String > mediumRuleNames;

  /// Field: specialTargets
  /// Each special makefile target is just a shorthand for a user-specified set of target files to make.
  map< String, vec< filename_t > > specialTargets;
  
  vec< filename_t > targetsToMake;

  /// Make command to use to run pipeline 
  String makeCommand;

  Bool dryRun;
  Bool lnForCp;
  Bool ignoreMetaData;
  Bool continueAfterError;
  Bool diagnostics;

  vec< filename_t > ignoreMissingMetadataFor;

  Bool ignoreInitialPartialOutputs;

  /**
     Field: ignoreChangesToTheseCommandArgs

     When deciding whether the command used to create a given
     file has changed ( and requires us to rebuild that file ),
     ignore the values of these command-line arguments.
     (They're simply erased from the command line before comparing
     the old and the new commands.)
   */
  vec< String > ignoreChangesToTheseCommandArgs;

  /**
     Field: dontUpdateTargetsOf

     Commands whose targets we will not update -- we will assume
     the targets already exist, and ignore the rule for this command.
   */
  set< shellcmd_t > dontUpdateTargetsOf;

  /**
     Field: forceTargetsOf

     Update targets of these commands even if they're up-to-date.
   */
  set< shellcmd_t > forceTargetsOf;

  vec< filename_t > forceTheseTargets;

  /**
     Field: restoreTargetsOf

     Commands whose targets we will restore from .incomplete files and
     then quit.
   */
  set< shellcmd_t > restoreTargetsOf;
  
  // private methods

  filename_t AbbrevFile( filename_t fn ) const;
  shellcmd_t AbbrevCommand( shellcmd_t command, String ignoreDef = "" ) const;
  filename_t CanonicalTarget( filename_t fn ) const;
  void MakePipelineDirectory(const String& directory) const; 

  void HandleTargetCaching( vec< filename_t > targetsToMake );
  
  static String EscapeForMake( String s );
  static String EscapeForShellEcho( String s );
  shellcmd_t CanonicalizeCommand( shellcmd_t cmd, Bool abbrevCmd = True ) const;


  void WriteDependencyGraph( filename_t dotFileName,
			     Bool hideTransitiveEdges ,
			     const vec< filename_t > *limitToTheseFiles,
			     const vec< rule_id_t > *limitToTheseRules,
			     const map< filename_t, url_t > *file2url,
			     const map< rule_id_t, url_t > *rule2url,
			     String graphName ) const;
  

  
  /**
     MethodDecl: RemovePartialOutputs

     Remove those files which might not be fully made.
     Any suspect file name F will be moved to F.incomplete.
  */
  void RemovePartialOutputs();

  /**
     MethodDecl: MakeInfo

     Given the full pathname of a target file,
     returns the basename of the corresponding metainfo
     files.  For each file somepath/F, the metainfo
     files are somepath/makeinfo/F.*,
     so here we return somepath/makeinfo/F
  */
  static filename_t MakeInfo( filename_t target );
  
  void ForceTargets( vec< filename_t > forceTargets ) const;
  void RemoveForcedTargets() const;
  Bool IsRuleRelevant( vec< MakeRule >::const_iterator r ) const;
   
};  // class MakeMgr



#endif
// #ifndef __INCLUDE_system_MiscUtil_h
