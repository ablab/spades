///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cctype>
#include <sstream>
#include "ParseSet.h"
#include "Set.h"
#include <list>
#include "String.h"
#include "system/Assert.h"
#include "system/MiscUtil.h"
#include "system/Assert.h"
#include "system/System.h"
#include "system/ParsedArgs.h"
#include "system/HTMLUtils.h"
#include "system/Utils.h"

static Bool makeOnlyIfNeeded = True;

void SetMakeOnlyIfNeeded( Bool newMakeOnlyIfNeeded ) {
  makeOnlyIfNeeded = newMakeOnlyIfNeeded;
}

/**
    Function: NeedToMake

    Given a set of target files and a set of source files,
    returns True if any of the target files is missing or
    or older than one of the source files.

    Halts with a fatal error if any of the source files
    are missing.
*/
Bool NeedToMake( const vec< filename_t >& targetFiles,
		 const vec< filename_t >& sourceFiles ) {
  ForceAssert( targetFiles.nonempty() );

  if ( !makeOnlyIfNeeded )
    return True;

  typedef int modtime_t;
  modtime_t latestSourceModTime = 0;

  for ( int i = 0; i < sourceFiles.isize(); i++ ) {
    if ( !IsRegularFile( sourceFiles[i] ) ) {
      cout << " Source file " << sourceFiles[i] << " missing!" << endl;
      cout << " Needed to make one of: " << endl << targetFiles << endl;
    }
    ForceAssert( IsRegularFile( sourceFiles[i] ) );
  }

  int latestSourceIdx = 0;
  for ( int i = 0; i < sourceFiles.isize(); i++ ) {
    modtime_t thisSourceModTime = LastModified( sourceFiles[i] );
    if ( i == 0  ||  thisSourceModTime > latestSourceModTime ) {
      latestSourceModTime = thisSourceModTime;
      latestSourceIdx = i;
    }
  }

  if ( sourceFiles.nonempty() )
    cout << " latest modified source is " << sourceFiles[ latestSourceIdx ] << endl;
  
  for ( int i = 0; i < targetFiles.isize(); i++ ) {
    if ( !IsRegularFile( targetFiles[i] ) ||
	 sourceFiles.nonempty()  &&  LastModified( targetFiles[i] ) < latestSourceModTime ) {
      cout << " Need to make " << targetFiles[i] << endl;
      if ( !IsRegularFile( targetFiles[i] ) )
	cout << " because it does not exist." << endl;
      else
	cout << " because its mod time is off." << endl;
      return True;
    }
  }

  cout << " The following files do not need to be made: " << targetFiles << endl;

  return False;
  
}

/**
   Function: Make

   Make target files from source files using the specified
   command.  The command is run only if needed,
   unless runOnlyIfNeeded==False.
*/
void Make( const vec< filename_t >& targetFiles,
	   const vec< filename_t >& sourceFiles,
	   const vec< shellcmd_t >& commands,
	   Bool runOnlyIfNeeded ) {
  if ( !runOnlyIfNeeded ||
       NeedToMake( targetFiles, sourceFiles ) ) {
    for ( int i = 0; i < commands.isize(); i++ ) {
      cout << " Running: " << commands[i] << endl;
      SystemSucceed( commands[i] );
      cout << " Finished: " << commands[i] << endl;
    }
  }
  ForceAssert( !NeedToMake( targetFiles, sourceFiles ) );
}

void Make( const vec< filename_t >& targetFiles,
	   const vec< filename_t >& sourceFiles,
	   shellcmd_t command,
	   Bool runOnlyIfNeeded ) {
  Make( targetFiles, sourceFiles, MkVec( command ), runOnlyIfNeeded );
}


/**
   Function: PrependDir

   Prepend a String to a vec of Strings -- for example, prepends a dirname to each filename.
*/
vec< filename_t > PrependDir( const filename_t& dir, const vec< filename_t >& files ) {
  vec< filename_t > result;
  for ( int i = 0; i < files.isize(); i++ ) {
    result.push_back( dir + "/" + files[i] );
  }
  return result;
}

/**
   FuncDecl: AbsToRelPaths

   Convert two absolute paths to a relative path.
*/

void AbsToRelPaths(const String& abs_source, const String& abs_dest, String& rel_source) {

  const String source_dir = abs_source.RevBefore("/");
  const String source_name = abs_source.RevAfter("/");

  unsigned int max_common_path_length = MIN(abs_source.size(), abs_dest.size());
  unsigned int pos = 0;
  while ( (pos < max_common_path_length) && (abs_dest[pos] == abs_source[pos]))
    ++pos;

  unsigned int common_end = abs_source.substr(0,pos).PosRev("/") + 1;

  String pathup = "";
  unsigned int up = abs_dest.substr(common_end).Freq("/");
  for (unsigned int i = 0; i < up; ++i)
    pathup += "../";

  String pathdown = "";
  if (common_end < source_dir.size())
    pathdown = source_dir.substr(common_end) + "/";

  rel_source = pathup + pathdown + source_name;
}


static vec< filename_t > RemoveEmpty( const vec< filename_t >& fnames ) {
  vec< filename_t > nonempty_fnames( fnames );
  EraseIf( nonempty_fnames, &String::empty );
  return nonempty_fnames;
}

/**
   Function: FilesIn

   Make a vector of filenames in the given dir.
 */
vec< filename_t > FilesIn( const filename_t& dir, const filename_t& fn1 ) {
  return PrependDir( dir, RemoveEmpty( MkVec( fn1 ) ) );
}
vec< filename_t > FilesIn( const filename_t& dir, const filename_t& fn1, const filename_t& fn2 ) {
  return PrependDir( dir, RemoveEmpty( MkVec( fn1, fn2 ) ) );
}
vec< filename_t > FilesIn( const filename_t& dir, const filename_t& fn1, const filename_t& fn2,
			   const filename_t& fn3) {
  return PrependDir( dir, RemoveEmpty( MkVec( fn1, fn2, fn3 ) ) );
}
vec< filename_t > FilesIn( const filename_t& dir, const filename_t& fn1, const filename_t& fn2,
			   const filename_t& fn3, const filename_t& fn4 ) {
  return PrependDir( dir, RemoveEmpty( MkVec( fn1, fn2, fn3, fn4 ) ) );
}
vec< filename_t > FilesIn( const filename_t& dir, const filename_t& fn1, const filename_t& fn2,
			   const filename_t& fn3, const filename_t& fn4, const filename_t& fn5 ) {
  return PrependDir( dir, RemoveEmpty( MkVec( fn1, fn2, fn3, fn4, fn5 ) ) );
}
vec< filename_t > FilesIn( const filename_t& dir, const filename_t& fn1, const filename_t& fn2,
			   const filename_t& fn3, const filename_t& fn4, const filename_t& fn5,
			   const filename_t& fn6 ) {
  return PrependDir( dir, RemoveEmpty( MkVec( fn1, fn2, fn3, fn4, fn5, fn6 ) ) );
}
vec< filename_t > FilesIn( const filename_t& dir, const filename_t& fn1, const filename_t& fn2,
			   const filename_t& fn3, const filename_t& fn4, const filename_t& fn5,
			   const filename_t& fn6, const filename_t& fn7 ) {
  return PrependDir( dir, RemoveEmpty( MkVec( fn1, fn2, fn3, fn4, fn5, fn6, fn7 ) ) );
}
vec< filename_t > FilesIn( const filename_t& dir, const filename_t& fn1, const filename_t& fn2,
			   const filename_t& fn3, const filename_t& fn4, const filename_t& fn5,
			   const filename_t& fn6, const filename_t& fn7, const filename_t& fn8 ) {
  return PrependDir( dir, RemoveEmpty( MkVec( fn1, fn2, fn3, fn4, fn5, fn6, fn7, fn8 ) ) );
}
vec< filename_t > FilesIn( const filename_t& dir, const filename_t& fn1, const filename_t& fn2,
			   const filename_t& fn3, const filename_t& fn4, const filename_t& fn5,
			   const filename_t& fn6, const filename_t& fn7, const filename_t& fn8,
			   const filename_t& fn9 ) {
  return PrependDir( dir, RemoveEmpty( MkVec( fn1, fn2, fn3, fn4, fn5, fn6, fn7, fn8, fn9 ) ) );
}


void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1 ) {
  ForceAssert( IsDirectory( dir1 ) && IsDirectory( dir2 ) );
  CpIfNewer( dir1 + "/" + fn1, dir2 + "/" + fn1 );
}

void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2 ) {
  ForceAssert( IsDirectory( dir1 ) && IsDirectory( dir2 ) );
  CpIfNewer( dir1 + "/" + fn1, dir2 + "/" + fn1 );
  CpIfNewer( dir1 + "/" + fn2, dir2 + "/" + fn2 );
}

void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2, filename_t fn3 ) {
  ForceAssert( IsDirectory( dir1 ) && IsDirectory( dir2 ) );
  CpIfNewer( dir1 + "/" + fn1, dir2 + "/" + fn1 );
  CpIfNewer( dir1 + "/" + fn2, dir2 + "/" + fn2 );
  CpIfNewer( dir1 + "/" + fn3, dir2 + "/" + fn3 );
}

void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2, filename_t fn3, filename_t fn4 ) {
  ForceAssert( IsDirectory( dir1 ) && IsDirectory( dir2 ) );
  CpIfNewer( dir1 + "/" + fn1, dir2 + "/" + fn1 );
  CpIfNewer( dir1 + "/" + fn2, dir2 + "/" + fn2 );
  CpIfNewer( dir1 + "/" + fn3, dir2 + "/" + fn3 );
  CpIfNewer( dir1 + "/" + fn4, dir2 + "/" + fn4 );
}

void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2, filename_t fn3, filename_t fn4, filename_t fn5 ) {
  ForceAssert( IsDirectory( dir1 ) && IsDirectory( dir2 ) );
  CpIfNewer( dir1 + "/" + fn1, dir2 + "/" + fn1 );
  CpIfNewer( dir1 + "/" + fn2, dir2 + "/" + fn2 );
  CpIfNewer( dir1 + "/" + fn3, dir2 + "/" + fn3 );
  CpIfNewer( dir1 + "/" + fn4, dir2 + "/" + fn4 );
  CpIfNewer( dir1 + "/" + fn5, dir2 + "/" + fn5 );
}

void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2, filename_t fn3, filename_t fn4, filename_t fn5, filename_t fn6 ) {
  ForceAssert( IsDirectory( dir1 ) && IsDirectory( dir2 ) );
  CpIfNewer( dir1 + "/" + fn1, dir2 + "/" + fn1 );
  CpIfNewer( dir1 + "/" + fn2, dir2 + "/" + fn2 );
  CpIfNewer( dir1 + "/" + fn3, dir2 + "/" + fn3 );
  CpIfNewer( dir1 + "/" + fn4, dir2 + "/" + fn4 );
  CpIfNewer( dir1 + "/" + fn5, dir2 + "/" + fn5 );
  CpIfNewer( dir1 + "/" + fn6, dir2 + "/" + fn6 );
}

void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2, filename_t fn3, filename_t fn4, filename_t fn5, filename_t fn6,
		filename_t fn7 ) {
  ForceAssert( IsDirectory( dir1 ) && IsDirectory( dir2 ) );
  CpIfNewer( dir1 + "/" + fn1, dir2 + "/" + fn1 );
  CpIfNewer( dir1 + "/" + fn2, dir2 + "/" + fn2 );
  CpIfNewer( dir1 + "/" + fn3, dir2 + "/" + fn3 );
  CpIfNewer( dir1 + "/" + fn4, dir2 + "/" + fn4 );
  CpIfNewer( dir1 + "/" + fn5, dir2 + "/" + fn5 );
  CpIfNewer( dir1 + "/" + fn6, dir2 + "/" + fn6 );
  CpIfNewer( dir1 + "/" + fn7, dir2 + "/" + fn7 );
}

void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2, filename_t fn3, filename_t fn4, filename_t fn5, filename_t fn6,
		filename_t fn7, filename_t fn8 ) {
  ForceAssert( IsDirectory( dir1 ) && IsDirectory( dir2 ) );
  CpIfNewer( dir1 + "/" + fn1, dir2 + "/" + fn1 );
  CpIfNewer( dir1 + "/" + fn2, dir2 + "/" + fn2 );
  CpIfNewer( dir1 + "/" + fn3, dir2 + "/" + fn3 );
  CpIfNewer( dir1 + "/" + fn4, dir2 + "/" + fn4 );
  CpIfNewer( dir1 + "/" + fn5, dir2 + "/" + fn5 );
  CpIfNewer( dir1 + "/" + fn6, dir2 + "/" + fn6 );
  CpIfNewer( dir1 + "/" + fn7, dir2 + "/" + fn7 );
  CpIfNewer( dir1 + "/" + fn8, dir2 + "/" + fn8 );
}

void CpIfNewer( dirname_t dir1, dirname_t dir2, filename_t fn1, filename_t fn2, filename_t fn3, filename_t fn4, filename_t fn5, filename_t fn6,
		filename_t fn7, filename_t fn8, filename_t fn9 ) {
  ForceAssert( IsDirectory( dir1 ) && IsDirectory( dir2 ) );
  CpIfNewer( dir1 + "/" + fn1, dir2 + "/" + fn1 );
  CpIfNewer( dir1 + "/" + fn2, dir2 + "/" + fn2 );
  CpIfNewer( dir1 + "/" + fn3, dir2 + "/" + fn3 );
  CpIfNewer( dir1 + "/" + fn4, dir2 + "/" + fn4 );
  CpIfNewer( dir1 + "/" + fn5, dir2 + "/" + fn5 );
  CpIfNewer( dir1 + "/" + fn6, dir2 + "/" + fn6 );
  CpIfNewer( dir1 + "/" + fn7, dir2 + "/" + fn7 );
  CpIfNewer( dir1 + "/" + fn8, dir2 + "/" + fn8 );
  CpIfNewer( dir1 + "/" + fn9, dir2 + "/" + fn9 );
}

MakeRule::MakeRule( const vec< filename_t >& targetFiles_, const vec< filename_t >& sourceFiles_,
		    shellcmd_t command_, filename_t saveOutputTo_, String comment_,
		    String shortRuleName_, String mediumRuleName_ ):
  targetFiles( targetFiles_ ), sourceFiles( sourceFiles_ ),
  command( command_ ), saveOutputTo( saveOutputTo_ ), comment( comment_ ),
  shortRuleName( shortRuleName_ ), mediumRuleName( mediumRuleName_ ) {
}


//
// Methods of class: MakeMgr
//

void MakeMgr::AddRule( vec< filename_t > targetFiles,
		       vec< filename_t > sourceFiles,
		       shellcmd_t command,
		       String comment,
		       filename_t saveOutputTo ) {
  if ( !preparedForRules )
    cout << "Please call MakeMgr::PrepareForRules() before "
      "adding any rules." << endl;
  ForceAssert( preparedForRules );

  ForceAssert( targetFiles.nonempty() );

  vec< filename_t > targetFilesExpanded = ExpandFileNames( ApplySubsts( targetFiles ) );
  for ( int i = 0; i < targetFilesExpanded.isize(); i++ ) {
    if ( STLContains( preExistingFiles, targetFilesExpanded[ i ] ) ) {
      StrongWarning( " Ignoring rule that makes " + AbbrevFile( targetFilesExpanded[ i ] )
		     + " because we were told to assume this file exists." );
      return;
    }
  }

  

  rule_id_t thisRuleId = rules.size();

  String shortRuleName;
  shellcmd_t cmdName = WhiteSpaceFree( command.SafeBefore(" ") );
  for ( int i = 0; i < cmdName.isize(); i++ )
    if ( isupper( cmdName[i] ) )
      shortRuleName += cmdName[i];
  if ( shortRuleName.empty() )
    shortRuleName = cmdName;

  if ( shortRuleName.empty() )
    shortRuleName = "rule_" + ToString( thisRuleId );

  String candidateShortRuleName = shortRuleName;
  int candidateShortRuleSuffix = 1;
  while ( STLContains( shortRuleNames, candidateShortRuleName ) )
    candidateShortRuleName = shortRuleName + "-" + ToString( ++candidateShortRuleSuffix );
  shortRuleNames.insert( shortRuleName = candidateShortRuleName );

  String mediumRuleName = cmdName;
  if ( mediumRuleName.empty() )
    mediumRuleName = shortRuleName;

  String candidateMediumRuleName = mediumRuleName;
  int candidateMediumRuleSuffix = 1;
  while ( STLContains( mediumRuleNames, candidateMediumRuleName ) )
    candidateMediumRuleName = mediumRuleName + "-" + ToString( ++candidateMediumRuleSuffix );
  ForceAssert( mediumRuleNames.find( candidateMediumRuleName ) == mediumRuleNames.end() );
  mediumRuleNames.insert( mediumRuleName = candidateMediumRuleName );

  if ( STLContains( dontUpdateTargetsOf, mediumRuleName ) ) {
    StrongWarning( " Treating targets of " + mediumRuleName + " as pre-existing." );
    preExistingFiles.insert( targetFilesExpanded.begin(), targetFilesExpanded.end() );
    return;
  }

  if ( getenv("DBG_MAKE") ) {
    cout << " Creating rule: " << endl;
    cout << " targets: " << ExpandFileNames( ApplySubsts( targetFiles ) ) << endl;
    cout << " sources: " << ExpandFileNames( ApplySubsts( sourceFiles ) ) << endl;
    cout << " command: " << ApplySubsts( command ) << endl;
  }

  rules.push( ExpandFileNames( ApplySubsts( targetFiles ) ),
	      ExpandFileNames( ApplySubsts( sourceFiles ) ),
	      ApplySubsts( command ),
	      ApplySubsts( saveOutputTo ),
	      comment,
	      shortRuleName,
	      mediumRuleName
	      );
  
  const vec< filename_t >& targets = rules.back().Targets();
  for ( int i = 0; i < targets.isize(); i++ ) {
    filename_t tf = targets[i];
    if ( STLContains( target2rule, tf ) ) {
      cout << " ERROR: conflicting rules for making target " << tf << endl;
      ForceAssert( False );
    } else {
      target2rule.insert( make_pair( tf, thisRuleId ) );
    }
  }
  
}  // MakeMgr::AddRule()

void MakeMgr::RemoveForcedTargets() const {
  ForceTargets( forceTheseTargets );
  For_( MakeRule, r, rules )
    if ( STLContains( forceTargetsOf, r->MediumRuleName() ) )
      ForceTargets( r->Targets() );
}


struct cmp_str_len: public binary_function<const String&, const String&, Bool> {
public:
  Bool operator() ( const String& s1, const String& s2 ) const { return s1.size() < s2.size(); }
};


/**
   Method: CheckTargetLocks

   For each target that we plan to make, check that it is not already being made.
   If it is, print an error message and return.
   If it isn't, create a lock 
 */

void MakeMgr::SetTargets( vec< filename_t > defaultTargets,
			  vec< filename_t > excludeTargets ) {
  
  defaultTargets = ExpandFileNames( ApplySubsts( defaultTargets ) );
  excludeTargets = ExpandFileNames( ApplySubsts( excludeTargets ) );
  vec< filename_t > defaultTargetsWithoutExcludes;
  set_difference( defaultTargets.begin(), defaultTargets.end(),
		  excludeTargets.begin(), excludeTargets.end(),
		  back_inserter( defaultTargetsWithoutExcludes ) );
  defaultTargets = defaultTargetsWithoutExcludes;
  targetsToMake.clear();
  vec< filename_t > active = defaultTargets;
  set< filename_t > seen(defaultTargets.begin(),defaultTargets.end());
  while ( active.nonempty() ) {
    filename_t current = active.back();
    active.pop_back();
    targetsToMake.push_back( current );
    map< filename_t, rule_id_t >::const_iterator it = target2rule.find( current );
    if ( it != target2rule.end() ) {
      vec< filename_t > deps = rules[ it->second ].Sources();
      for (int i = 0; i < deps.isize(); i++) {
	if (seen.insert(deps[i]).second)
	  active.push_back( deps[i]);
      }
    }
  }
  UniqueSort( targetsToMake );
  PRINT( targetsToMake );
}

void MakeMgr::MakePipelineDirectory(const String& directory) const {
  if ( !IsDirectory( directory ) ) {
    cout << " Creating " << AbbrevFile( directory ) << endl;
    if (dryRun) {
      cout << " (not actually creating directory because this is a dry run)" << endl;
    } else {
      SystemSucceed( "mkdir -p " +  directory );
      ForceAssert( IsDirectory(  directory ) );
    }
  }
}


/**
   Method: TestRules
   
   Tests dependency and target list for each rule.

*/
int MakeMgr::TestRules(String top_dir, String rule, String mode, bool erase ) {

  cout << Date() << " : Testing Rules" << endl << endl;

  if ( targetsToMake.empty() ) {
    cout << " *** NO TARGETS TO MAKE!!" << endl;
    return 1;
  }

  // Determine directories
  while (top_dir.Contains("//"))
    top_dir.GlobalReplaceBy( "//", "/" );
  if (top_dir.EndsWith("/"))
    top_dir = top_dir.SafeBefore("/");
  const String base_dir = top_dir.SafeBeforeLast("/");  
  const String orig_name = top_dir.SafeAfterLast("/");
  const String orig_dir = base_dir + "/" + orig_name;

  bool global_error = false;

  bool single_rule = (rule != "ALL" && rule != "all" && rule != "All");

  // Test rules, one by one
  for ( size_t ruleIdx = 0; ruleIdx < rules.size(); ruleIdx++ ) {

    // Track problems
    bool error = false;
    bool warning = false;

    // Grab rule and skip if it is a copy or link rule
    const MakeRule& r = rules[ ruleIdx ];
    String rule_name = r.MediumRuleName();
    String cmd_name = r.CommandName();
    if ( cmd_name == "cp" || cmd_name == "ln" )  // skip copy or link rules
      continue;

    // Test all rules or just one?
    if (single_rule) {
      if (rule.IsInt() && ( rule.Int() != static_cast<int>(ruleIdx)))  // test rule N only
	continue;
      else if (!rule.IsInt() && (rule != rule_name)) // test rule COMMAND only
	continue;
    }

    cout << "Rule: " << ruleIdx << ",  Name: " << rule_name << endl;
    cout << "============================================================================"  <<  endl;

    // Skip rules that might not have the dependencies to run (unless explicitly asked for)
    bool skip = true;
    for (size_t i = 0; i < r.Targets().size() && skip; i++)
      skip = !Member(targetsToMake, r.Targets()[i]);
    if (skip) {
      if (single_rule) 
	cout << " Warning! - rule not used in current pipeline" << endl << endl;
      else {
	cout << " Skipping test - rule not used in current pipeline" << endl << endl << endl;
	continue;
      }
    }

    String test_name = orig_name + "/test_rules/" + ToString(ruleIdx) + "." + rule_name;
    String test_dir = base_dir + "/" + test_name;
    String log_file = test_dir + "/test_rule.log";

    if (mode != "list") {
      // Create testing directory for this rule
      if (IsDirectory(test_dir)) {
	cout << "ERROR, existing test directory found: " << test_dir << endl;
	return 1;
      }
      Mkpath(test_dir);
    }
    
    // Translate dependencies to test directory
    vec< filename_t > sources = r.Sources();
    for (size_t i = 0; i < sources.size(); i++) 
      sources[i].ReplaceBy(orig_name, test_name);

    // Translate targets to test directory
    vec< filename_t > targets = r.Targets();
    for (size_t i = 0; i < targets.size(); i++)
      targets[i].ReplaceBy(orig_name, test_name);

    // Translate output target (if it exists) to test directory
    String out_target = (r.SaveOutputTo() == "" ? "" :
			 r.SaveOutputTo().ReplaceBy(orig_name, test_name));

    // Translate command to test directory
    String command =  r.Command();
    command.GlobalReplaceBy(orig_name, test_name);

    if (mode != "list") {
      // Symlink dependencies to test directory
      for (size_t i = 0; i < sources.size(); i++) {
	Mkpath(Dirname(sources[i]));
	if (IsRegularFile(r.Sources()[i]))
	  Symlink(r.Sources()[i], sources[i]);
	else {
	  cout << "ERROR - unable to find original file:" << endl << r.Sources()[i] << endl;
	  return 1;
	}	
      }
    }

    if (mode == "prepare" || mode == "list")
      cout << " Command:" << endl << command << endl;

    if (mode == "prepare") {
      cout << " Test directory prepared:" << endl << test_dir << endl;
      continue;
    }

    // Run command
    if (mode != "list") {
      int status = System( command + (single_rule ? "" : " > " + log_file + " 2>&1") );
      if (single_rule) cout << endl;
      if (status != 0) {
	if (!single_rule) 
	  cout << " Command failed, possible missing dependency - see test_rule.log file." << endl;
	error = true;
      } else {
	
	if (mode == "module")
	  continue;

	// Check all targets were made (skipping output target)
	vec <String> missing;
	for (size_t i = 0; i < targets.size(); i++) {
	  if (( targets[i] != out_target) && !IsRegularFile(targets[i]) )
	    missing.push_back(targets[i]);
	}
	if (!missing.empty()) {
	  error = true;
	  cout << " Missing targets:" << endl;
	  for (size_t i = 0; i < missing.size(); i++)
	    cout << "  " << missing[i].SafeAfter(test_name + "/") << endl;
	}
	
	// Remove targets and dependencies
	for (size_t i = 0; i < targets.size(); i++) 
	  if (IsRegularFile(targets[i]) ) 
	    Remove( targets[i]);
	for (size_t i = 0; i < sources.size(); i++) 
	  Remove( sources[i]);
	
	// Check for undeclared targets 
	vec <string> undeclared;
	vec <String> files = AllFiles(test_dir);
	list <String> to_check(files.begin(), files.end());
	while (to_check.size() != 0) {
	  String& item = to_check.front();
	  if (IsDirectory(test_dir + "/" + item) ) {
	    files = AllFiles(test_dir + "/" + item);
	    for (size_t i = 0; i < files.size(); i++)
	      to_check.push_back(item + "/" + files[i]);
	  } else if  (item != "test_rule.log")
	    undeclared.push_back(item);
	  to_check.pop_front();	
	}
	if (!undeclared.empty()) {
	  warning = true;
	  cout << " Undeclared targets:" << endl;
	  for (size_t i = 0; i < undeclared.size(); i++) 
	    cout << "  " << undeclared[i] << endl;
	}

	// Check for unused dependencies (very time consuming)
	if (mode == "full") {
	  vec <String> unused;
	  for (size_t i = 0; i < sources.size(); i++) {
	    // clear and refresh test directory
	    System("rm -r " + test_dir);
	    for (size_t link_index = 0; link_index < sources.size(); link_index++) {
	      Mkpath(Dirname(sources[link_index]));
	    if (i != link_index)
	      Symlink(r.Sources()[link_index], sources[link_index]);
	    }
	    if ( System( command + " > " + log_file + " 2>&1") == 0 ) {
	      // Check all targets were made
	      bool missing_target = false;
	      for (size_t target_index = 0; target_index < targets.size(); target_index++) 
		missing_target |=  !IsRegularFile(targets[target_index]);
	      if (!missing_target)
		unused.push_back(sources[i]);
	    }
	  }
	  if (!unused.empty()) {
	    warning = true;
	    cout << " Unused dependencies:" << endl;
	    for (size_t i = 0; i < unused.size(); i++)
	      cout << "  " << unused[i].SafeAfter(test_name + "/") << endl;
	  }
	}
      }
    }
    // Rule testing complete
    
    // On error or warning, report expected dependencies and targets
    if (error || warning || (mode == "list") ) {
      cout << " Targets: " << endl;
      for (size_t i = 0; i < targets.size(); i++) 
	cout << "  " << targets[i].SafeAfter(test_name + "/") << endl;
      cout << " Dependencies: " << endl;
      for (size_t i = 0; i < sources.size(); i++) 
	cout << "  " << sources[i].SafeAfter(test_name + "/") << endl;
    }
    if (error) 
      cout << " FAILED!" << endl;
    else if (warning)  
      cout << " Warning!" << endl;
    else
      cout << " Passed!" << endl;
    
    // Remove test directory
    if (erase && mode != "list")
      System("rm -r " + test_dir);
    
    cout << endl << endl;
    global_error |= error;
  }

  if (global_error) {
    cout << "At least one rule failed - see above for more details." << endl;
    return 1;
  }

  return 0;
}



/**
   Method: RunMake
   
   Run the make process.

   Implementation:h

   Generates a make makefile, then runs make.

*/
int MakeMgr::RunMake( 
		      int MAX_PARALLEL ) {
  //  char makeFileName [128];
  //  strcpy( makeFileName, "makeXXXXXX" );

  cout << Date() << " : Preparing makefile." << endl << endl;

  RemoveForcedTargets();


  if ( targetsToMake.empty() ) {
    cout << " *** NO TARGETS TO MAKE!!" << endl;
    return 1;
  }

  //
  // In each dir that contains one or more of our targets,
  // make sure there is a makeinfo/ subdir.
  //

  vec< dirname_t > targetDirs;
  for ( map< filename_t, rule_id_t >::const_iterator it = target2rule.begin();
	it != target2rule.end(); it++ )
    targetDirs.push_back ( Dirname( it->first ) );
  UniqueSort( targetDirs );
  for ( int i = 0; i < targetDirs.isize(); i++ ) {
    dirname_t makeInfoDir = targetDirs[i] + "/makeinfo";
    MakePipelineDirectory( makeInfoDir );
  }
  
  // Create logging directory
  MakePipelineDirectory(logDirectory);
  
  const String makeLog(logDirectory + "/" + "make.log");
  
  // Summary pipeline log file
  const String summaryLog(logDirectory + "/" + "summary.log");
  const String pipeToLog(" >> " + summaryLog);

  temp_file makeFileName( "/tmp/MakeMgr_XXXXXX" );

  {
    Ofstream( mf, makeFileName );

    //
    // Define makefile variables for commonly-used directories,
    // as specified by earlier calls to AddDirAbbrev().
    // This makes the makefile, and the dependency graph,
    // more human-readable.
    //
    if ( !dirAbbrevs.empty() ) {
      mf << "\n\n";

      vec< dirname_t > abbreviatedDirs = Reverse( abbreviatedDirsRevByLen );

      for ( int i = 0; i < abbreviatedDirs.isize(); i++ ) {
	dirname_t d = abbreviatedDirs[i];
	dir2abbrev_t::const_iterator it = dirAbbrevs.find( d );

	dirname_t abbreviatedDef = AbbrevCommand( d, it->second );
	mf << it->second << "=" <<
	   ( abbreviatedDef.Contains( "$(" + it->second + ")" ) ? it->first : abbreviatedDef )
	   << "\n";
      }
	
      mf << "\n\n";
    }

    
    if ( targetsToMake.nonempty() ) {
      
      vec< filename_t > canonTargets;
      for ( int i = 0; i < targetsToMake.isize(); i++ )
	canonTargets.push_back( AbbrevFile( CanonicalTarget( targetsToMake[i] ) ) );
      UniqueSort( canonTargets );
      mf << "default: ";
      for ( int i = 0; i < canonTargets.isize(); i++ )
	mf << " " << canonTargets[i];
      mf << "\n\n";
    }

    mf << ".PHONY: all clean\n";
    mf << "\nall: ";
    for ( int i = 0; i < rules.isize(); i++ )
      mf << AbbrevFile( rules[i].PrimaryTarget() ) << " ";
    mf << "\n\n";

    
    mf << "\nclean: \n\trm -f ";
    for ( int i = 0; i < rules.isize(); i++ ) {
      for ( int j = 0; j < rules[i].Targets().isize(); j++ )
	mf << AbbrevFile( rules[i].Targets()[j] ) << " ";
    }
    mf << "\n\n";
    
    for ( int ruleIdx = 0; ruleIdx < rules.isize(); ruleIdx++ ) {
      const MakeRule& r = rules[ ruleIdx ];
      const vec< filename_t >& sources = r.Sources();
      const vec< filename_t >& targets = r.Targets();
      const filename_t primaryTarget = r.PrimaryTarget();
      ForceAssert( targets.nonempty() );
      const String RS = "[" + r.ShortRuleName() + "] ";
      const String echo = EscapeForShellEcho( "\t@echo " + RS );
      const String echodate = EscapeForShellEcho("\t@echo ") + RS + " `date` : ";
      mf << "\n";
      if ( targets.isize() > 1 ) {
	mf << "# ";
	for ( int i = 1; i < targets.isize(); i++ )
	  mf << AbbrevFile( targets[i] ) << ", ";
	mf << endl;
      }
      mf << AbbrevFile( primaryTarget ) << ": ";
      vec< filename_t > depsList;
      for ( int sourceIdx = 0; sourceIdx < sources.isize(); sourceIdx++ )
	depsList.push_back( AbbrevFile( CanonicalTarget( sources[ sourceIdx ] ) ) );
      UniqueSort( depsList );
      for ( int i = 0; i < depsList.isize(); i++ )
	mf << depsList[i] << " ";
      mf << "\n";

      shellcmd_t cmdName = r.CommandName();
      
      String outputSavedTo(" ");
      
      if ( r.SaveOutputTo().nonempty() )
	outputSavedTo += AbbrevFile( r.SaveOutputTo() ) + " ";
//       for ( vec< filename_t >::const_iterator target = r.Targets().begin();
// 	    target != r.Targets().end();  target++ )
// 	outputSavedTo += AbbrevFile( MakeInfo( *target ) + ".out." + cmdName ) + " ";
      outputSavedTo += AbbrevFile( MakeInfo( *(r.Targets().begin()) ) + ".out." + cmdName ) + " ";
      outputSavedTo += AbbrevFile( logDirectory + "/" + r.MediumRuleName() + ".out" ) + " ";


      String cmdSavedTo(" ");
      for ( vec< filename_t >::const_iterator target = r.Targets().begin();
	    target != r.Targets().end();  target++ ) {
	cmdSavedTo += AbbrevFile( MakeInfo( *target ) + ".cmd" ) + " ";
      }

      mf << "\t@rm -f ";
      for ( int i = 0; i < targets.isize(); i++ )
	mf << AbbrevFile( targets[i] ) << " ";
      mf << outputSavedTo << " " << cmdSavedTo;
      mf << "\n";

      if ( r.SaveOutputTo().nonempty() ) {
	mf << echo << "Saving output of " << cmdName << " to " << EscapeForShellEcho( AbbrevFile( r.SaveOutputTo() ) ) << "\n";
      }

      if ( cmdName != "cp" && cmdName != "ln" ) {
	//
	// Tell the user what we're doing.
	//
	mf << echodate << "Starting command " << cmdName << "\n";
	mf << echo << "\n";
	mf << echo << " === " << EscapeForShellEcho( r.Comment() ) << " ===\n";
	mf << echo << "\n";
	mf << echo << "Calling " << cmdName << " to create " << targets.isize() << " file\\(s\\):\n";
	mf << echo << "\n";
	for ( int i = 0; i < targets.isize(); i++ )
	  mf << echo << EscapeForShellEcho( AbbrevFile( targets[i] ) )<< endl;
	mf << echo << "\n";
	mf << echo << "from " << sources.isize() << " file\\(s\\):\n";
	mf << echo << "\n";
	for ( int i = 0; i < sources.isize(); i++ )
	  mf << echo << EscapeForShellEcho( AbbrevFile( sources[i] ) )<< endl;
	mf << echo << "\n";
      }

      mf << echo << EscapeForShellEcho( AbbrevCommand( r.Command() ) ) << endl;
      mf << echodate << r.MediumRuleName() << " Starting" << pipeToLog << endl;

      /////////////////
      // Run the command itself!  We have so much other boilerplate here,
      // the command we're trying to run to update the targets
      // seems to have gotten lost...
      ////////////////
      mf << "\t@" << AbbrevCommand( r.Command() );
      ///////////////
      
      if ( cmdName != "cp" && cmdName != "gzip" && cmdName != "ln" ) {
        if ( diagnostics ) {
          mf << " MM=True _MM_INTERVAL=10 _MM_SUMMARY=False _MM_OUT=" 
	    + AbbrevFile( MakeInfo( primaryTarget ) ) + ".mm." + cmdName + " ";
        }

        mf << " TEE=\"" << outputSavedTo << "\"  2>&1 \n";
      } else {
	mf << "\n";
      }

      if ( cmdName != "cp" && cmdName != "gzip" && cmdName != "ln" ) 
	mf << echodate << " Validating targets" << "\n";
      
      // Check that all targets got made.  If not, immediately stop.
      for ( int i = 0; i < targets.isize(); i++ )
	mf << "\t@if [ ! -f " << AbbrevFile( targets[i] ) << " ]; then echo Program " 
	   << cmdName << " failed to make target " 
	   <<  EscapeForShellEcho( AbbrevFile( targets[i] ) ) << "; "
	   << "@echo `date` : " << r.MediumRuleName() << " failed to make target "
	   << EscapeForShellEcho( AbbrevFile( targets[i] ) ) << pipeToLog << "; "
	   << "exit 1; fi\n";

      // touch all targets to ensure they all have the same modification date
      if ( cmdName != "cp" && cmdName != "gzip" && cmdName != "ln" ) {	
	mf << "\t@touch -c ";
	for ( int i = 0; i < targets.isize(); i++ )
	  mf << AbbrevFile( targets[i] ) << "  "; 
	mf << "\n";
      }
	
      mf << echodate << r.MediumRuleName() << " Finished" << pipeToLog << endl;

      //
      // For each target F, write out the command used to create the target,
      // and the dependencies, into makeinfo/F.cmd
      //
      // We first write this out for the primary target,
      // then copy to cmd files for the other targets.
      //
      filename_t cmdFile = AbbrevFile( MakeInfo( primaryTarget ) ) + ".cmd";
      String toCmdFile = " >> " + cmdFile + "\n";
      mf << "\t@echo " << EscapeForShellEcho( AbbrevCommand( r.Command() ) ) << toCmdFile;
      mf << "\t@echo \\[sources\\]" << toCmdFile;
      for ( int i = 0; i < sources.isize(); i++ )
	mf << "\t@echo " << EscapeForShellEcho( AbbrevFile( sources[i] ) ) << toCmdFile;
      mf << "\t@echo \\[targets\\]" << toCmdFile;
      for ( int i = 0; i < targets.isize(); i++ )
	mf << "\t@echo " << EscapeForShellEcho( AbbrevFile( targets[i] ) ) << toCmdFile;
      mf << "\t@echo \\[comment\\]" << toCmdFile;
      mf << "\t@echo " << ( r.Comment().empty() ? String( "none" ) : EscapeForShellEcho( r.Comment() ) ) << toCmdFile;
      mf << "\t@echo \\[end\\]" << toCmdFile;
      for ( int i = 0; i < targets.isize(); i++ )
	if ( targets[i] != primaryTarget )
	  mf << "\t@cp " << cmdFile << " " << AbbrevFile( MakeInfo( targets[i] ) ) << ".cmd\n";
      mf << echodate << "Finished command " << cmdName << "\n";
      
      if ( primaryTarget != r.SaveOutputTo() ) {
	mf << "\n.DELETE_ON_ERROR: " << AbbrevFile( primaryTarget ) << "\n";
      }
    }  // for each rule
    
    
    {
      
      //
      // Check the dependency graph for errors
      //
      
      // Identify makefile targets which we do not know how to make.
      // These files must alread exist for us to be able to complete
      // the make.
      
      
      vec< pair< filename_t, rule_id_t > > targetsWithNoRule;
      for ( int i = 0; i < rules.isize(); i++ ) {
	for ( int j = 0; j < rules[i].Sources().isize(); j++ ) {
	  filename_t src = rules[i].Sources()[j];
	  if ( !STLContains( target2rule, src ) &&
	       !STLContains( preExistingFiles, src ) )
	    targetsWithNoRule.push( src, i );
	}
      }
      UniqueSort( targetsWithNoRule );
      
      if ( targetsWithNoRule.nonempty() )
	cout << "\n\n targets with no rule: " << targetsWithNoRule.size() << endl;
      for ( int i = 0; i < targetsWithNoRule.isize(); i++ ) {
	cout << AbbrevFile( targetsWithNoRule[i].first ) <<
	  " ( needed by " + rules[ targetsWithNoRule[i].second ].MediumRuleName() +
	  " to make " << AbbrevFile( rules[ targetsWithNoRule[i].second ].PrimaryTarget() )
	     << " when " << rules[ targetsWithNoRule[i].second ].Comment()
	     << " )" << endl;
      }      
      cout << "\n" << endl;
      ForceAssert( targetsWithNoRule.empty() );
      
      for( set< filename_t >::const_iterator it = preExistingFiles.begin();
	   it != preExistingFiles.end(); it++ ) {
	if ( !IsRegularFile( *it )  &&
	     BinMember( targetsToMake, *it ) ) {
	  cout << " Error: file " << *it << " is supposed to already exist, but doesn't." << endl;
	  ForceAssert( IsRegularFile( *it ) );
	}
      }

      // Clean up incomplete builds: for any rule that did not complete and
      // build _all_ its targets, delete the targets to force the rule
      // to be re-run.  This most often happens when we change a program
      // to output a file that it did not output before.
    }
  }


  if ( !ignoreInitialPartialOutputs ) 
    RemovePartialOutputs();

  if (outFilename != "") {
    System( "cp " + makeFileName + " " + outFilename );
    cout << "Copy of makefile in: " << outFilename << endl;
  }

  String command = makeCommand + (dryRun ? " -n " : " ") + ToString(" -r ") +
    (continueAfterError ? " -k " : "") + ToString(" -f ") + makeFileName +
    (MAX_PARALLEL > 1 ? " -j" + ToString( MAX_PARALLEL ) : ToString("") ) ;


  // Run make, but first test to see if targets are already up to date.
  int retCode = System( command + " -q"  );

  if (retCode == 0) {
    // Targets are all up to date (1 - nothing to do)
    return 1;
  } else if (retCode == 1) {
    // Some targets are not up to date, really need to make...
    cout << endl << Date() << " : make process starting." << endl << endl;
    //    int retCode =  System( command + (logDirectory == "" ? "" : " 2>&1 | tee " + makeLog) );
    //    int retCode =  System( "mkfifo pipe ; tee " + makeLog + " < pipe &  " + command + " > pipe");
    retCode =  System( command );
    cout << endl << Date() << " : make process finished." << endl;
    RemovePartialOutputs();
    return (retCode == 0 ? 0 : 2);  // 0 = success;  2 = error occured
  }

  if (retCode == 2) {
    // Problem occured during make, write message to summary file
    OfstreamMode(summary_log_file, summaryLog, ios::app );
    summary_log_file << "\n  *** Pipeline encountered an error, "
		     <<"see full logs for detailed error messages. ***\n" << endl;
  }
  
  // Return exit code (0 - success, 2 - error, (and previously: 1 - nothing to do) )
  return retCode;
 
}

void MakeMgr::AddCpRule( filename_t target, filename_t source ) {
  if ( lnForCp ) {
    String dest_dir, dest_name, rel_source;
    AbsToRelPaths(source, target, rel_source);
    AddRule( MkVec( target ), MkVec( source ),
	      "ln -sf " + rel_source + " " + target );
  } else AddRule( MkVec( target ), MkVec( source ),
		"cp " + source + " " + target );
}

void MakeMgr::AddCpRule( dirname_t targetDir, dirname_t sourceDir, const vec< filename_t >& fnames ) {
  vec< filename_t > fnamesExpanded= ExpandFileNames( ApplySubsts( fnames ) );
  for ( int i = 0; i < fnamesExpanded.isize(); i++ )
    AddCpRule( targetDir + "/" + fnamesExpanded[i], sourceDir + "/" + fnamesExpanded[i] );
}

void MakeMgr::AddCpRule( dirname_t targetDir, filename_t targetStart, 
			 dirname_t sourceDir, filename_t sourceStart,
			 const filename_t  commonEnd ) {
  vec< filename_t > commonEndExpanded = ExpandFileNames( MkVec ( ApplySubsts( commonEnd ) ) );
  for ( unsigned int i = 0; i < commonEndExpanded.size(); ++i )
    AddCpRule( targetDir + "/" + targetStart + commonEndExpanded[i], 
	       sourceDir + "/" + sourceStart + commonEndExpanded[i] );
}


void MakeMgr::AddCpRule( dirname_t targetDir, dirname_t sourceDir, filename_t fn1 ) {
  AddCpRule( targetDir, sourceDir, MkVec( fn1 ) );
}


filename_t MakeMgr::ApplySubsts( filename_t fn ) const {
  for ( map< String, String >::const_iterator it = substs.begin(); it != substs.end(); it++ )
    fn.GlobalReplaceBy( it->first, it->second );
  return fn;
  
}

vec< filename_t > MakeMgr::ApplySubsts( const vec< filename_t >& fnames ) const {
  vec< filename_t > result;
  for ( int i = 0; i < fnames.isize(); i++ )
    result.push_back( ApplySubsts( fnames[i] ) );
  return result;
}


filename_t MakeMgr::AbbrevFile( filename_t fn ) const {
  //dir2abbrev_t::const_iterator abbrev = dirAbbrevs.find( Dirname( fn ) );
  //String result = abbrev != dirAbbrevs.end() ? "$(" + abbrev->second + ")/" + Basename(fn) : fn;
  return ApplySubsts( AbbrevCommand( fn ) );
}

void MakeMgr::PrepareForRules() {
  ForceAssert( !preparedForRules );
  
  for ( dir2abbrev_t::const_iterator it = dirAbbrevs.begin();
	it != dirAbbrevs.end(); it++ )
    abbreviatedDirsRevByLen.push_back( it->first );
  ReverseSort( abbreviatedDirsRevByLen,
	       cmp_str_len() );

  

  preparedForRules = True;
}

shellcmd_t MakeMgr::AbbrevCommand( shellcmd_t cmd, String ignoreDef ) const {
  for ( int i = 0; i < abbreviatedDirsRevByLen.isize(); i++ ) {
    dirname_t d = abbreviatedDirsRevByLen[i];
    dir2abbrev_t::const_iterator it = dirAbbrevs.find( d );
    if ( it != dirAbbrevs.end() && it->second != ignoreDef )
      cmd.GlobalReplaceBy( d,
			   "$(" + it->second + ")" );
  }

  return cmd;
}

filename_t MakeMgr::CanonicalTarget( filename_t fn ) const {
  map< filename_t, rule_id_t >::const_iterator it =
    target2rule.find( fn );
  filename_t result = it == target2rule.end() ? fn : rules[ it->second ].PrimaryTarget();
  return result;
}


void MakeMgr::WriteDependencyGraph( filename_t dotFileName,
				    Bool hideTransitiveEdges,
				    const vec< filename_t > *limitToTheseFiles,
				    const vec< rule_id_t > *limitToTheseRules,
				    const map< filename_t, String > *file2url,
				    const map< rule_id_t, String > *rule2url,
				    String graphName ) const {

  // Construct a digraph showing the connections between the files and the commands.

  // make a list of all files
  vec< filename_t > allFiles;

  for ( int i = 0; i < rules.isize(); i++ ) {
    allFiles.append( rules[i].Sources() );
    allFiles.append( rules[i].Targets() );
  }
  UniqueSort( allFiles );
  
  // assign a node name for each file
  map< filename_t, int > fname2id;
  map< filename_t, String > fname2nodeName;
  for ( int i = 0; i < allFiles.isize(); i++ ) {
    fname2id[ allFiles[i] ] = i;
    fname2nodeName[ allFiles[i] ] = "file_" + BaseAlpha( i );
  }

  matrix<Bool> A( rules.size() + allFiles.size(), rules.size() + allFiles.size(), False );

  for ( int ruleId = 0; ruleId < rules.isize(); ruleId++ ) {
    const MakeRule& r = rules[ ruleId ];
    if ( !limitToTheseRules || BinMember( *limitToTheseRules, ruleId ) ) {

      for ( int i = 0; i < r.Sources().isize(); i++ ) {
	if ( !limitToTheseFiles || BinMember( *limitToTheseFiles, r.Sources()[i] ) )
	  A( rules.size() + fname2id[ r.Sources()[i] ], ruleId ) = True;
      }
      for ( int i = 0; i < r.Targets().isize(); i++ ) {
	if ( !limitToTheseFiles || BinMember( *limitToTheseFiles, r.Targets()[i] ) )
	  A( ruleId, rules.size() + fname2id[ r.Targets()[i] ] ) = True;
      }
    }
  }
  
  digraph G( A );

  Ofstream( df, dotFileName );
  df << "digraph " << graphName << " {\n";
  df << "    rankdir=\"LR\";\n";
  
  int edgesSkipped = 0;

  vec< filename_t > filesUsed;
    
  for ( rule_id_t ruleIdx = 0; ruleIdx < rules.isize(); ruleIdx++ ) {
    if ( limitToTheseRules && !BinMember( *limitToTheseRules, ruleIdx ) )
      continue;
    
      const MakeRule& r = rules[ ruleIdx ];
    
      String ruleNodeName = "rule_" + BaseAlpha( ruleIdx );
      shellcmd_t cmdToolTip;
      String theCmd = r.Comment().nonempty() ? r.Comment() : 
	AbbrevCommand( r.Command() );
      for ( int i = 0; i < theCmd.isize(); i++ ) {
	char c = theCmd[i];
	if ( c == '\"' )
	  cmdToolTip += '\\';
	cmdToolTip += c;
      }
      
      df << ruleNodeName << " [ label=\"" << r.MediumRuleName() << "\", shape=parallelogram, fillcolor=green, tooltip=\"" +
	cmdToolTip + "\" ";
      map< rule_id_t, String >::const_iterator it;
      if ( rule2url && ( it = rule2url->find( ruleIdx ) ) != rule2url->end() )
	df << ", URL=\"" << it->second << "\"";
      df << " ];\n";
      
      for ( int i = 0; i < r.Sources().isize(); i++ ) {

	if ( limitToTheseFiles && !BinMember( *limitToTheseFiles, r.Sources()[i] ) )
	  continue;
	
	Bool transitivePathFound = False;
	
	if ( hideTransitiveEdges ) {
	  
	  // if there is a longer path, don't show this edge
	  vec< vec< vrtx_t > > paths;
	  vrtx_t v = rules.size() + fname2id[ r.Sources()[ i ] ];
	  vrtx_t w = ruleIdx;
	  if ( A( v, w ) && G.From( v ).nonempty() && G.To( w ).nonempty() ) {

	    if ( !G.From( v ).solo() && !G.To( w ).solo() ) {
	      if ( G.AllPaths( v, w, paths ) ) {
		for ( int k = 0; k < paths.isize() && !transitivePathFound; k++ ) {
		  const vec< vrtx_t >& path = paths[k];
		  ForceAssertEq( path.front(), v );
		  ForceAssertEq( path.back(), w );
		  if ( path.size() > 2 )
		    transitivePathFound = True;
		}
	      }
	    }
	    
	  }
	}
	
	if ( !transitivePathFound ) {
	  if ( !limitToTheseFiles  ||  BinMember( *limitToTheseFiles, r.Sources()[i] ) ) {
	    df << fname2nodeName[ r.Sources()[i] ] << " -> " <<
	      ruleNodeName << " [ tailtooltip=\"" << AbbrevFile( r.Sources()[i] ) << "\", tailURL=\"http://www.mit.edu\" ];\n";
	    filesUsed.push_back( r.Sources()[i] );
	  }
	} else {
	  //cout << " skipping edge: " << endl;
	  //cout << fname2nodeName[ r.Sources()[i] ] << " -> " <<
	  //ruleNodeName << " [ tailtooltip=\"" << AbbrevFile( r.Sources()[i] ) << "\", tailURL=\"http://www.mit.edu\" ];\n";
	  edgesSkipped++;
	}
      }  // for all sources
      for ( int i = 0; i < r.Targets().isize(); i++ ) {
	if ( !limitToTheseFiles || BinMember( *limitToTheseFiles, r.Targets()[i] ) ) {
	  df << ruleNodeName << " -> " << fname2nodeName[ r.Targets()[i] ] << " [ tailtooltip=\"" << r.MediumRuleName() << "\", tailURL=\"http://www.boston.com\" ];\n";
	  filesUsed.push_back( r.Targets()[i] );
	}
      }
    }  // for each rule


  UniqueSort( filesUsed );
  for ( int i = 0; i < allFiles.isize(); i++ )
    if ( BinMember( filesUsed, allFiles[i] ) &&
	 ( !limitToTheseFiles || BinMember( *limitToTheseFiles, allFiles[i] ) ) ) {
      df << fname2nodeName[ allFiles[i] ] << " [ label=\"" << AbbrevFile( allFiles[i] ) << "\", shape=oval";
      map< filename_t, String >::const_iterator it;
      if ( file2url && (it = file2url->find( allFiles[i] )) != file2url->end() )
	df << ", URL=\"" << it->second << "\"";
      map< filename_t, rule_id_t >::const_iterator it2 =
	target2rule.find( allFiles[i] );

      while ( it2 != target2rule.end()  &&  rules[ it2->second ].CommandName() == "cp" &&
	      rules[ it2->second ].Sources().solo() )
	it2 = target2rule.find( rules[ it2->second ].Sources().front() );
      
      df << ", tooltip=\"";
      if ( it2 == target2rule.end() )
	df << "assumed to exist";
      else {
	const MakeRule& fileCreatorRule = rules[ it2->second ];

	df << fileCreatorRule.MediumRuleName() << ": ";
	df << ( fileCreatorRule.Comment().nonempty() ? fileCreatorRule.Comment() :
		AbbrevCommand( fileCreatorRule.Command() ) );
      }
  
  
      df << "\" ];\n";
    }
  
  df << "}\n";
}


/**
   Method: HandleTargetCaching
   
   For each target that we need to build and that already exists,
   check that we're building it with the same command as before.  If not,
   move the current target to the cache, and if the cache has a version
   of this target that's built with the currently requested parameters,
   then retrieve that target from the cache.
   
   Implementation notes:
   
   Now, the target we retrieve may end up getting overwritten because
   some of its prereqs have changed.  In that case we should perhaps
   just delete it from the cache.  For the moment, let make decide.
   
   Now, about the order of arguments: perhaps not the whole command-line
   matters for the result, and the argument order most likely does not
   matter.  So, if it's one of Arachne commands, we should perhaps
   call ParsedArgs to parse it and then write down the canonical
   version.  We can add this later.   For now, the most common
   idiom will be where we copy one of two versions of a file
   to a common "join point" name.
*/
void MakeMgr::HandleTargetCaching( vec< filename_t > targetsToMake ) {
#if 0   
  // Identify the list of targets we'll be making.
  
  if ( targetsToMake.empty() ) return;
  
  vec< filename_t > allTargetsToMake( targetsToMake );
  vec< filename_t > active;
  active.push( targetsToMake.front() );
  while ( active.nonempty() ) {
    filename_t current = active.back();
    active.pop_back();
    allTargetsToMake.push_back( current );
    map< filename_t, rule_id_t >::const_iterator it =
      target2rule.find( fn );
    if ( it != target2rule.end() ) {
      active.append( rules[ it->second ].Sources() );
    }
    
    UniqueSort( allTargetsToMake );
    
    for ( int i = 0; i < alTargetsToMake; i++ ) {
      filename_t target = allTargetsToMake[i];
      if ( IsRegularFile( target ) ) {
	
      }
    }
  }
#endif  
}
  
void MakeMgr::ForceTargets( vec< filename_t > forceTargets ) const {
  forceTargets = ExpandFileNames( ApplySubsts( forceTargets ) );
  PRINT( forceTargets );
  for ( int i = 0; i < forceTargets.isize(); i++ ) {
    if ( forceTargets[i] != "all" && forceTargets[i] != "clean" ) {
      filename_t f = CanonicalTarget( forceTargets[i] ) ;
      if ( IsRegularFile( f ) ) {
	dirname_t makeinfoDir = Dirname( MakeInfo( f ) );
	if ( !IsDirectory( makeinfoDir ) )
	  SystemSucceed( "mkdir -p " + makeinfoDir );
	filename_t f_moved = MakeInfo(f) + ".incomplete";
	cout << " Moving " << AbbrevFile( f ) << " to "
	     << AbbrevFile( f_moved ) << " to force it to be remade." << endl;
	SystemSucceed( "mv -f " +  f + " " + f_moved );
      }
    }
  }
}

void MakeMgr::AddSpecialTarget( String specialTarget, vec< filename_t > filesInTarget ) {
  ForceAssertNe( specialTarget, "all" );
  ForceAssertNe( specialTarget, "clean" );
  specialTargets.insert( make_pair( specialTarget, filesInTarget ) );
}


String MakeMgr::EscapeForMake( String s ) {
  String r;
  for ( int i = 0; i < s.isize(); i++ ) {
    if ( s[i] == '$' )
      r += '$';
    r += s[i];
  }
  return r;
}

shellcmd_t MakeMgr::CanonicalizeCommand( shellcmd_t cmd, Bool abbrevCmd ) const {
  shellcmd_t r = cmd;

  // remove any of the ignored args.
  for ( int i = 0; i < ignoreChangesToTheseCommandArgs.isize(); i++ ) {
    String ignoredArg = ignoreChangesToTheseCommandArgs[i];

    int pos;
    while ( ( pos = r.Position( " " + ignoredArg + "=" ) ) != -1 ) {
      int nextSpacePos = pos + 1;
      while ( nextSpacePos < r.isize()  &&  !isspace( r[ nextSpacePos ] ) )
	nextSpacePos++;
      r.erase( pos+1, nextSpacePos - pos );
    }
  }
  
  DeleteLeadingWhiteSpace( r );
  DeleteTrailingWhiteSpace( r );

  // compress multiple adjacent whitespace chars into single whitespace char.
  // also, convert any whitespace chars to a plain space.
  String r_canonSpaces;
  for ( int i = 0; i < r.isize(); i++ ) {
    char c = r[i];
    if ( !( i > 0 && isspace( r[i-1] ) && isspace( c ) ) )
      r_canonSpaces += ( isspace( c ) ? ' ' : c );
  }
  r = r_canonSpaces;
  
  return abbrevCmd ? AbbrevCommand( r ) : r;
}



String MakeMgr::EscapeForShellEcho( String s ) {
  String r;
  static const char *charsToEscape = "\"\'\\{}()[]$`;";
  for ( int i = 0; i < s.isize(); i++ ) {
    if ( strchr( charsToEscape, s[i] ) )
      r += '\\';
    r += s[i];
  }
  return EscapeForMake( r );
}

void MakeMgr::RestoreTargets() {
  //  UniqueSort( restoreTargetsOf );
  int numRestored = 0;
  For_( MakeRule, r, rules ) {
    if ( STLContains( restoreTargetsOf, r->MediumRuleName() ) ) {
      For_( filename_t, target, r->Targets() ) {
	filename_t mvFrom = *target, mvTo = MakeInfo( *target ) + ".incomplete";
	if ( IsRegularFile( mvTo ) ) {
	  cout << " Moving " << AbbrevFile( mvTo ) << " back to " << AbbrevFile( mvFrom ) << endl;
	  SystemSucceed( "mv -f " +  mvTo + " " + mvFrom );
	  numRestored++;
	}
      }
    }
  }
  cout << " Restored " << numRestored << " targets, exiting..." << endl;
  exit( 0 );

}

Bool MakeMgr::IsRuleRelevant( vec< MakeRule >::const_iterator r ) const {
  For_( filename_t, target, r->Targets() ) {
    if ( BinMember( targetsToMake, *target ) ) {
      //cout << "Relevant Rule: " << r->MediumRuleName() << endl;
      return True;
    }
  }

  //cout << "Irrelevant Rule: " << r->MediumRuleName() << endl;
  return False;
}

/**
   Method: RemovePartialOutputs
   
   Remove those files which might not be fully made.
   Any suspect file name D/F will be moved to D/makeinfo/F.incomplete.
*/
void MakeMgr::RemovePartialOutputs() {

  UniqueSort( ignoreMissingMetadataFor );
  
  For_ ( MakeRule, r, rules ) {

    /**

       Check that the rule ran to completion.
       If that's not the case, remove ALL targets of the rule.

       When a rule has run successfully to completion, the following is true:

       - for each target F, there is a makeinfo/F.cmd file create AFTER the
       target, giving the command used to create the target, as well as
       the list of sources used by the command, and the list of targets
       written by it.


       The rule should be remade if
    
       - any source has been updated since the rule was run
       ( make will check this )
	 
    */

    // First, check that the rule makes some target
    // that we plan to make.

    if ( !IsRuleRelevant( r ) )
      continue;
    
    int numTargetsExisting = 0;
    vec< filename_t > missing, incomplete, commandChanged;
    vec< filename_t > missingButNotNeeded;
    vec< pair< shellcmd_t, shellcmd_t > > commandChangedFromTo;
    vec< pair< filename_t, filename_t > > newSourcesFound;
    For_ ( filename_t, target, r->Targets() ) {
      if ( IsRegularFile( *target ) ) {
	numTargetsExisting++;
	filename_t cmdFileName = MakeInfo( *target ) + ".cmd";
	if ( ignoreMetaData ) {
	  // don't check metadata file - assume all is well
	} else if ( !IsRegularFile( cmdFileName ) ) {
	  if ( BinMember( ignoreMissingMetadataFor, *target ) ) {
	    // cout << " Ignoring missing metadata for " << AbbrevFile( *target ) << endl;
	  } else
	    incomplete.push_back( *target );
	} else if ( IsOlder( cmdFileName, *target ) ) {
	  cout << "Command metadata file is older than target " << AbbrevFile(*target) << endl;
	  incomplete.push_back( *target );
	} else {
	  Ifstream( cmdFile, cmdFileName );
	  shellcmd_t oldCmd;
	  getline( cmdFile, oldCmd );
	  String canonOld = CanonicalizeCommand( oldCmd );
	  String canonNew = CanonicalizeCommand( r->Command() ) ;
	  if ( canonOld != canonNew &&
	       // also try canonicalizing without abbreviations, to deal with legacy files
	       CanonicalizeCommand( oldCmd, False ) != CanonicalizeCommand( r->Command(), False ) ) {

#if 0	    
	    cout << " canonicalize of " << endl;
	    cout << oldCmd << endl;
	    cout << " is " << endl;
	    cout << canonOld << endl;
	    cout << " canonicalize of " << endl;
	    cout << r->Command() << endl;
	    cout << " is " << endl;
	    cout << canonNew << endl;

	    PRINT2( canonOld.size(), canonNew.size() );

	    for ( int i = 0; i < canonOld.isize() && i < canonNew.isize(); i++ ) {
	      if ( canonOld[i] != canonNew[i] ) {
		PRINT5( i, char(canonOld[i]), int(canonOld[i]), char(canonNew[i]), int(canonNew[i]) );
	      }
	    }
	    if ( canonNew.size() > canonOld.size() ) {
	      for ( int i = canonOld.isize(); i < canonNew.isize(); i++ )
		PRINT3( i, char( canonNew[i] ), int( canonNew[i] ) );
	    }
#endif	    
	      
	    commandChanged.push_back( *target );
	    commandChangedFromTo.push( canonOld, canonNew );
	  }
	  if ( cmdFile.good() ) {
	    String sectionHeading;
	    getline( cmdFile, sectionHeading );
	    if ( sectionHeading == "[sources]" ) {
	      vec< filename_t > sourcesBefore;
	      while( cmdFile.good() ) {
		String s;
		getline( cmdFile, s );
		if ( s.empty() || s.StartsWith( "[" ) )
		  break;
		sourcesBefore.push_back( s );
	      }
	      UniqueSort( sourcesBefore );
	      for ( int i = 0; i < r->Sources().isize(); i++ ) {
		filename_t sourceNow = AbbrevFile( r->Sources()[i] );
		if ( !BinMember( sourcesBefore, sourceNow ) &&
		     // also try the unabbreviated version, to deal
		     // with legacy files
		     !BinMember( sourcesBefore, r->Sources()[i] ) )
		  newSourcesFound.push( *target, sourceNow );
	      }
	    }
	  }
	}
      } else  { // !IsRegularFile( *target )
	if ( !BinMember( targetsToMake, *target ) ) {
	  missingButNotNeeded.push_back( *target );
	} else 
	  missing.push_back( *target );
      }
    }
    if ( numTargetsExisting > 0 ) {
      for ( int k = 0; k < missingButNotNeeded.isize(); k++ ) {
	// cout << " File " << AbbrevFile( missingButNotNeeded[k] ) << " is missing, but we do not need it." << endl;
      }
    }
    if ( numTargetsExisting > 0  &&
	 ( missing.nonempty()  ||  incomplete.nonempty() ||
	   commandChanged.nonempty() || newSourcesFound.nonempty() ) ) {
      if ( missing.nonempty() || incomplete.nonempty() ) {
	cout << "\n\n Looks like the command " << r->CommandName() << " did not complete: \n";
	if ( missing.nonempty() ) {
	  cout << " The following files are missing: " << missing.size() << endl;
	  for ( int i = 0; i < missing.isize(); i++ )
	    cout << AbbrevFile( missing[i] ) << endl;
	}
	if ( incomplete.nonempty() ) {
	  cout << " The following files are incomplete: " << incomplete.size() << endl;
	  for ( int i = 0; i < incomplete.isize(); i++ )
	    cout << AbbrevFile( incomplete[i] ) << endl;	
	}
      }
      if ( commandChanged.nonempty() ) {
	cout << "\n\nLooks like the command line of command " << r->MediumRuleName() << " has changed: \n";
	set< pair< shellcmd_t, shellcmd_t > > seen;
	for ( int i = 0; i < commandChanged.isize(); i++ ) {
	  if( !STLContains( seen, commandChangedFromTo[i] ) ) {
	    cout << " Command line to create " << AbbrevFile( commandChanged[i] ) << " changed from\n";
	    cout << commandChangedFromTo[i].first << endl;
	    cout << " to \n";
	    cout << commandChangedFromTo[i].second << endl;
	    cout << "\n";
	    seen.insert( commandChangedFromTo[i] );
	  }
	}
      }
      if ( newSourcesFound.nonempty() ) {
	cout << "\n\nLooks like the dependencies (source files) of the command " << r->MediumRuleName() << " have changed: \n";
	for ( int i = 0; i < newSourcesFound.isize(); i++ ) {
	  cout << "File \n" << newSourcesFound[i].first << "\nhas new dependency on\n" << newSourcesFound[i].second << "\n\n";
	}
      }

      cout << " Cleaning up all targets of this command..." << endl;
      for ( vec< filename_t >::const_iterator target = r->Targets().begin();
	    target != r->Targets().end(); target++ ) {
	filename_t mvFrom = *target, mvTo = MakeInfo( *target ) + ".incomplete";
	if ( IsRegularFile( mvFrom ) ) {
	  cout << " Moving " << AbbrevFile( mvFrom ) << " to " << AbbrevFile( mvTo ) << endl;
	  if ( dryRun )
	    cout << " (not actually moving because this is a dry run)" << endl;
	  else
	    SystemSucceed( "mv -f " +  mvFrom + " " + mvTo );
	}
      }
    }
  }
}


void MakeMgr::AddIgnoreChangesToCommandArgs( String commandArgs ) {
    ForceAssert( !preparedForRules );

    vec< String > commandArgSet;
    ParseStringSet( commandArgs, commandArgSet, true /* recurse */ );
    ignoreChangesToTheseCommandArgs.append( commandArgSet );
}


/**
   Method: MakeInfo
   
   Given the full pathname of a target file,
   returns the basename of the corresponding metainfo
   files.  For each file somepath/F, the metainfo
   files are somepath/makeinfo/F.*,
   so here we return somepath/makeinfo/F
*/
filename_t MakeMgr::MakeInfo( filename_t target ) {
  return Dirname( target ) + "/makeinfo/" + Basename( target );
}


/**
    Method: WritePipelineDescrAsHTML

    Create a webpage that describes the entire pipeline.
*/
void MakeMgr::WritePipelineDescrAsHTML( dirname_t dir ) const {
  if ( !IsDirectory( dir ) )
      SystemSucceed( "mkdir -p " + dir );
  ForceAssert( IsDirectory( dir ) );

  String cp_str = (lnForCp ? "ln" : "cp");

  Ofstream( f, dir + "/ruleDetail.html" );
  f << HTMLHead( "Pipeline" ) << endl;

  vec< filename_t > allFiles;


  for ( int i = 0; i < rules.isize(); i++ ) {
    allFiles.append( rules[i].Sources() );
    allFiles.append( rules[i].Targets() );
  }
  UniqueSort( allFiles );

  if ( parsed_args::ARGS ) {
    f << "<pre>\n";
    ((parsed_args *)parsed_args::ARGS)->PrintTheCommandPretty( f );
    f << "</pre>\n";
  
    f << "<hr>\n";
  }
  
  f << "<p>" << rules.size() << " rules to make " << target2rule.size() << " targets.</p>\n" << endl;
  f << "<p><a name=\"assumedToExist\"></a>" << preExistingFiles.size() << " files assumed to exist:</p>\n" << endl;
  f << "<ul>\n";
  for ( set< filename_t >::const_iterator it = preExistingFiles.begin(); it != preExistingFiles.end(); it++ )
    f << "<li>" << AbbrevFile( *it ) << "</li>" << endl;
  f << "</ul>\n";

  f << "<hr>\n";

  // for each target, make a list of rules that use that target as a source
  map< filename_t, vec< rule_id_t > > target2users;

  for ( rule_id_t i = 0; i < rules.isize(); i++ ) {
    for ( int j = 0; j < rules[i].Sources().isize(); j++ )
      target2users[ rules[i].Sources()[j] ].push_back( i );
  }
  

  {
    Ofstream( df, dir + "/ruleGraph.dot" );
    df << "digraph ruleGraph {\n";

    matrix<Bool> A( rules.size(), rules.size(), False );

    rule_id_t i = 0;
    For_( MakeRule, r, rules ) {

      if ( r->CommandName() != cp_str && IsRuleRelevant( r ) ) {
	//if ( !IsRuleRelevant( r ) ) continue;

	vec< rule_id_t > allSuccessors;

	vec< filename_t > targetsToUse( r->Targets() );
	while( targetsToUse.nonempty() ) {
	  filename_t targ = targetsToUse.back();
	  targetsToUse.pop_back();
	  map< filename_t, vec< rule_id_t > >::const_iterator usersIt = target2users.find( targ );
	  if ( usersIt != target2users.end() ) {
	    const vec< rule_id_t >& rids = usersIt->second;
	    for ( vec< rule_id_t >::const_iterator ridsIt = rids.begin(); ridsIt != rids.end(); ridsIt++ ) {
	      rule_id_t rid = *ridsIt;
	      if ( rules[ rid ].CommandName() == cp_str ) {
		targetsToUse.append( rules[ rid ].Targets() );
	      } else {
		if ( IsRuleRelevant( rules.begin() + rid ) )
		  allSuccessors.push_back( rid );
	      }
	    }
	  }
	}
      
	UniqueSort( allSuccessors );

	for ( int j = 0; j < allSuccessors.isize(); j++ )
	  A( i, allSuccessors[j] ) = True;

      }
      
      i++;
    } // for each rule

    digraph G( A );
    i = 0;
    For_ ( MakeRule, r, rules ) {
      if ( r->CommandName() != cp_str && IsRuleRelevant( r ) ) {
	//if ( !IsRuleRelevant( r ) ) continue;
	String ruleNodeName = "rule_" + BaseAlpha( i );
	shellcmd_t cmdToolTip;
	String theCmd = r->Comment().nonempty() ? r->Comment() : 
	  AbbrevCommand( r->Command() );
	for ( int k = 0; k < theCmd.isize(); k++ ) {
	  char c = theCmd[k];
	  if ( c == '\"' )
	    cmdToolTip += '\\';
	  cmdToolTip += c;
	}

	for ( int j = 0; j < G.From( i ).isize(); j++ ) {

	  vec< vec< vrtx_t > > paths;
	  Bool transitivePathFound = False;
	  if ( G.AllPaths( i, G.From(i)[j], paths ) ) {
	    for ( int k = 0; k < paths.isize() && !transitivePathFound; k++ ) {
	      const vec< vrtx_t >& path = paths[k];
	      if ( path.size() > 2 )
		transitivePathFound = True;
	    }
	  }

	  if ( !transitivePathFound )
	    df << "rule_" << BaseAlpha( i ) << " -> rule_" << BaseAlpha( G.From(i)[j] ) << " ;\n";

	}

	df << ruleNodeName << " [ label=\"" << r->MediumRuleName() << "\", shape=parallelogram, fillcolor=green, tooltip=\"" <<
	  cmdToolTip << "\", target=\"ruleDetail\", URL=\"ruleDetail.html#rule" << ToString(i) << "\" ];\n\n";

      }  // if the rule is not a cp rule

      i++;
    }  // for each rule
    
    df << "}\n";
  }

  {
    filename_t dotFileBase = dir + "/ruleGraph";
    shellcmd_t dotCmd = "dot -Tcmapx -o" + dotFileBase + ".map -Tgif -o" + dotFileBase + ".gif " + dotFileBase + ".dot";
    cout << " running " << dotCmd << endl;
    SystemSucceed( dotCmd );

    {
      Ofstream( ff, dir + "/index.html" );
      ff << HTMLHeadFrames( " Pipeline " );
      ff << "<frameset rows=\"75%,25%\" >" << endl;
      ff << "<frame name=\"ruleGraph\" id=\"ruleGraph\" src=\"ruleGraph.html\" >\n";
      ff << "<frame name=\"ruleDetail\" id=\"ruleDetail\" src=\"ruleDetail.html\">"<< endl;
      ff << "</frameset>"<< endl;
      ff << HTMLTailFrames();
    }


    {
      Ofstream( ff, dir + "/ruleGraph.html" );
    
      ff << "<img src=\"ruleGraph.gif\" USEMAP=\"#ruleGraph\" >\n";
      Ifstream (is, dotFileBase + ".map") ;
      String theMap = Slurp( is );
      ff << theMap << endl;
    }
  }



  map< filename_t, String > file2url;
  for ( int i = 0; i < allFiles.isize(); i++ ) {
    map< filename_t, rule_id_t >::const_iterator it2 = target2rule.find( allFiles[i] );
    while ( it2 != target2rule.end()  &&  rules[ it2->second ].CommandName() == cp_str &&
	    rules[ it2->second ].Sources().solo() )
      it2 = target2rule.find( rules[ it2->second ].Sources().front() );
    
    file2url.insert( make_pair( allFiles[i], it2 != target2rule.end() ? "#rule" + ToString( it2->second ) :
				ToString( "#assumedToExist" ) ) );
  }
  
  map< rule_id_t, String > rule2url;
  for ( rule_id_t i = 0; i < rules.isize(); i++ )
    rule2url.insert( make_pair( i, "#rule" + ToString( i ) ) );

  rule_id_t ruleNum = 0;
  for ( vec< MakeRule >::const_iterator r = rules.begin(); r != rules.end(); r++, ruleNum++ ) {
    f << "<hr>\n";

    f << "<a name=\"rule" << ruleNum << "\"></a>\n";

    f << "<p><b>" << r->MediumRuleName() << "</b>: " << r->Comment() << "</p>\n";
    
    vec< filename_t > relevantFiles;
    relevantFiles.append( r->Sources() );
    relevantFiles.append( r->Targets() );
    vec< rule_id_t > relevantRules;
    relevantRules.push_back( ruleNum );

    vec< filename_t > targetsToUse( r->Targets() );
    while( targetsToUse.nonempty() ) {
      filename_t targ = targetsToUse.back();
      targetsToUse.pop_back();
      map< filename_t, vec< rule_id_t > >::const_iterator usersIt = target2users.find( targ );
      if ( usersIt != target2users.end() ) {
	const vec< rule_id_t >& rids = usersIt->second;
	relevantRules.append( rids );
	for ( vec< rule_id_t >::const_iterator ridsIt = rids.begin(); ridsIt != rids.end(); ridsIt++ ) {
	  rule_id_t rid = *ridsIt;
	  if ( rules[ rid ].CommandName() == cp_str ) {
	    relevantFiles.append( rules[ rid ].Targets() );
	    targetsToUse.append( rules[ rid ].Targets() );
	  }
	}
      }
    }
    
    UniqueSort( relevantFiles );
    UniqueSort( relevantRules );

    filename_t dotFileBase = dir + "/rule" + ToString( ruleNum );
    WriteDependencyGraph( dotFileBase + ".dot", True /* hide transitive edges */,
			  &relevantFiles, &relevantRules, &file2url, &rule2url, "rulemap" + ToString( ruleNum) );
    shellcmd_t dotCmd = "dot -Tcmapx -o" + dotFileBase + ".map -Tgif -o" + dotFileBase + ".gif " + dotFileBase + ".dot";
    cout << " running " << dotCmd << endl;
    SystemSucceed( dotCmd );
    f << "<img src=\"rule" << ruleNum << ".gif\" USEMAP=\"#rulemap" << ruleNum << "\" >\n";
    Ifstream (is, dotFileBase + ".map") ;
    String theMap = Slurp( is );
    f << theMap << endl;

    f << "<p><code>" << AbbrevCommand( r->Command() ) << "</code></p>" << endl;
    
  }  // for each rule

  f << HTMLTail() << endl;
}  // WritePipelineDescrAsHTML




vec< filename_t> MakeMgr::ExpandFileNames( const vec< filename_t >& fnames ) const {
  vec< filename_t > result;
  for ( int i = 0; i < fnames.isize(); i++ ) {
    vec< filename_t > expanded;
    if ( fnames[i].empty() ) continue;
    ParseStringSet( fnames[i], expanded, true /* recurse */ );
    for ( int j = 0; j < expanded.isize(); j++ ) {

      if ( expanded[j] == "all" ) {
	for ( int k = 0; k < rules.isize(); k++ )
	  result.append( rules[k].Targets() );
      } else {
	map< String, vec< filename_t > >::const_iterator it =
	  specialTargets.find( expanded[j] );
	if ( it != specialTargets.end() )
	  result.append( ExpandFileNames( it->second ) );
	else
	  result.push_back( expanded[j] );
      }
      
    }
  }
  UniqueSort( result );
  return result;
}

