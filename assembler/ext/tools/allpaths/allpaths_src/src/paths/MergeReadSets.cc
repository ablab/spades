///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


/* MergeReadSets
 *
 * Merges two or more fastb files into one, also merging any associated qualb,
 * paths, pairing, and qltout files together too.
 *
 * If the kmer space differs between the various input paths, you should set
 * REPATH=True to re-path.
 *
 * If the input qltout files use different reference sequences, the output
 * qltout file will be nonsensical.
 *
 * If one of these file types does not exist in the input datasets, it will not be
 * created in the output dataset.
 *
 * Input files take the form: DIR/READS_IN.*
 * Output files, if created, will take the form: DIR/READS_OUT.*.
 *
 * Files merged are:
 *     fastb (required)
 *     qualb
 *     paths.k<K>  (only if K is given - also set REPATH if necessary)
 *     pairs/pairto
 *     qltout
 *
 *****************************************************************************/


#include "Basevector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "lookup/LookAlign.h"
#include "paths/ReadsToPathsCoreX.h"
#include "system/System.h" 
#include "util/ReadTracker.h"
#include "feudal/BinaryStream.h"
#include "math/IntDistribution.h"


static inline 
String Tag(String S = "MRS") { return Date() + " (" + S + "): "; } 



bool CheckFileSetExists(const vec<String>& heads, const String& suffix, bool verbose = false) {
  bool found_all = true;
  for (size_t i = 0; i < heads.size(); ++i)
    if (IsRegularFile(heads[i] + suffix) == false) {
      if (verbose) 
	cout << "Could not find file: " << heads[i] + suffix << endl;
      found_all = false;
    }
  return found_all;
}


// Merge the feudal files at <heads_in>.suffix into a single file at
// <head_out>.suffix.
void MergeFeudal(const String& head_out, const vec<String>& heads_in, const String& suffix, const vec<uint64_t> & sizes_check) {
  vec<String> full_in(heads_in.size());
  for (size_t i = 0; i < heads_in.size(); ++i) {
    full_in[i] = heads_in[i] + suffix;
    ForceAssertEq(MastervecFileObjectCount(full_in[i]), sizes_check[i]);
  }
  Remove(head_out + suffix);
  MergeMastervecs(head_out + suffix, full_in);
}



// Merge the alignment files at <heads_in>.qltout into a single file at
// <head_out>.qltout.  The dataset sizes are important, because qltout files
// do not know how many reads they contain.
void MergeQltout(const String& head_out, const vec<String>& heads_in, const vec<uint64_t> & sizes) {
  
  size_t offset = 0;
  Ofstream(out, head_out + ".qltout");
  
  // Loop over all of the input files.
  for (size_t i = 0 ; i < heads_in.size(); i++) {
    if (IsRegularFile(heads_in[i] + ".mapq"))
         Cp(heads_in[i] + ".mapq", head_out + ".mapq", True);
    
    // Load LookAligns from this input file.
    vec<look_align> aligns;
    LoadLookAligns(heads_in[i] + ".qltout", aligns);
    
    // Apply an offset to each of these LookAligns' read IDs, to account for the
    // merging.  Then write to the output file.
    for (size_t j = 0; j < aligns.size(); j++) {
      ForceAssertLt((uint64_t)aligns[j].query_id, sizes[i]);
      aligns[j].query_id += offset;
      aligns[j].PrintParseable(out);
      aligns[j].PrintReadableBrief(out);
    }
    
    offset += sizes[i];

  }
}


int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;
  CommandDoc("Merges two or more fastb files into one, but also merges associated qualb, "
	      "paths, pairing, and qltout files together too. Optionally recomputes paths if the kmer "
	      "space differs between the various input paths.");
  CommandArgument_StringSet_Doc(READS_IN,
    "Comma separated list of fastb files to merge.");
  CommandArgument_String_Doc(READS_OUT,
    "Merged fastb filename.");
  CommandArgument_String_OrDefault_Doc(DIR, "",
    "Location of <READS_IN> files to merge and the output files.");
  CommandArgument_Int_OrDefault_Doc(K, 0,
    "K size if merging paths or REPATHing");
  CommandArgument_Bool_OrDefault_Doc(REPATH, False,
    "Recompute paths instead of merging. Set to true if the existing paths don't share "
    "the same kmer space");
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use when REPATHing (use all available processors if set to 0)");
  CommandArgument_Bool_OrDefault_Doc(TRACK_READS, True,
    "Use ReadTracker to track source of merged reads");
  CommandArgument_Bool_OrDefault_Doc(MERGE_QUALS, True,
    "If available, merge qualb files.");
  CommandArgument_Bool_OrDefault_Doc(MERGE_PATHS, True,
    "If available, merge paths files.");
  CommandArgument_Bool_OrDefault_Doc(MERGE_PAIRS, True,
    "If available, merge pairs or pairto files.");
  CommandArgument_Bool_OrDefault_Doc(MERGE_ALIGNS, True,
    "If available, merge qltout files.  If the input qltout files use different reference sequences, the output qltout file will be nonsensical.");
  CommandArgument_Bool_OrDefault_Doc(MERGE_DISTRIBS, True,
    "If available, merge paired read separations distribution files.");
  CommandArgument_Bool_OrDefault_Doc(FORCE_MERGE_DISTRIBS, True,
    "merge paired read separations distribution files. Make them up if none available.");
  EndCommandArguments;
  
  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);

  String head_out = (DIR == "" ? "" : DIR + "/") + READS_OUT.SafeBefore(".fastb");
 
  String paths_suffix = ".paths.k" + ToString(K);
  
  // Expand input filenames, stripping off .fastb extension if necessary
 
  vec<String> heads_in(READS_IN.size());
  for (size_t i = 0; i < READS_IN.size(); ++i) {
    heads_in[i] = (DIR == "" ? "" : DIR + "/") + READS_IN[i].SafeBefore(".fastb");
    if (heads_in[i] == head_out)
      FatalErr("Input and output filenames conflict: " + READS_IN[i]);
  }
  const size_t n_sets = heads_in.size();


  // Check to see which files are available to merge.
  
  if (CheckFileSetExists(heads_in, ".fastb", true) == false)
    FatalErr("Fastb file missing - see above for details");
  
  bool found_paths    = MERGE_PATHS && CheckFileSetExists(heads_in, paths_suffix);
  bool found_quals    = MERGE_QUALS && CheckFileSetExists(heads_in, ".qualb");
  bool found_pairs    = MERGE_PAIRS && CheckFileSetExists(heads_in, ".pairs");
  bool found_pairto   = MERGE_PAIRS && CheckFileSetExists(heads_in, ".pairto");
  bool found_qltout   = MERGE_ALIGNS && CheckFileSetExists(heads_in, ".qltout");
  bool found_distribs = MERGE_DISTRIBS && CheckFileSetExists(heads_in, ".distribs");

  // Find the sizes of the input datasets.
  vec<size_t> sizes;
  for (size_t i = 0; i < READS_IN.size(); i++)
    sizes.push_back(MastervecFileObjectCount(heads_in[i] + ".fastb"));
  size_t size_sum = BigSum(sizes);
  
  cout << Tag() << "Merging Details:" << endl;
  if (DIR != "")
    cout << "In directory:" << endl 
	 << "  " << DIR << endl;
  cout << "Merging the following readsets:" << endl;
  for (size_t i = 0; i < READS_IN.size(); ++i) 
    cout << "  " << READS_IN[i].SafeBefore(".fastb") << ".*\t\t[n reads = "
	 << sizes[i] << "]" << endl;
  cout << "To create output readset:" << endl;
  cout << "  " << READS_OUT.SafeBefore(".fastb") << ".*\t\t\t[will have n reads = "
       << size_sum << "]" << endl;

  // ---- Merge fastb files
  
  cout << Tag() << "Merging fastb files" << endl;
  MergeFeudal(head_out, heads_in, ".fastb", sizes);
  
  // ---- Merge qualb files

  if (found_quals) {
    cout << Tag() << "Merging qualb files" << endl;
    MergeFeudal(head_out, heads_in, ".qualb", sizes);
  }
 



  // ---- Merge pairing info

  if (found_pairs || found_pairto || found_distribs) {

    vec< vec<String> > inLibNames(n_sets);
    PairsManager pairs;
    for (size_t i = 0; i < n_sets; ++i) {
      PairsManager pairs_loc;
      if (found_pairs) {
	pairs_loc.Read(heads_in[i] + ".pairs");
      } else {
	size_t n_reads = MastervecFileObjectCount(heads_in[i] + ".fastb");
	pairs_loc.ReadFromPairtoFile(heads_in[i] + ".pairto", n_reads);
      }
      ForceAssertEq(pairs_loc.nReads(), sizes[i]);
      pairs.Append(pairs_loc);
      inLibNames[i] = pairs_loc.getLibraryNames();
    }
  
    if (found_pairs || found_pairto) {  
      if (found_pairs) {
	cout << Tag() << "Writing merged pairs files" << endl;
	pairs.Write(head_out + ".pairs");
      }
      else {
	cout << Tag() << "Writing merged pairto files" << endl;
	vec<read_pairing> pairings = pairs.convert_to_read_pairings();
	WritePairs(pairings, pairs.nReads(), head_out, false);
      }
    }




    
    // ---- Merge separation distributions

    if (found_distribs || FORCE_MERGE_DISTRIBS) {
      cout << Tag() << "Merging read separation distribution files" << endl;
      size_t nlibs = pairs.nLibraries();
      vec<IntDistribution> distribs(nlibs);

      for (size_t i = 0 ; i < n_sets; i++) {   // Loop over all of read sets.

        const size_t n_libs = inLibNames[i].size();

	if (IsRegularFile(heads_in[i] + ".distribs")) {

	  vec<IntDistribution> distribs_loc;
	  const String fileIn = heads_in[i] + ".distribs";
	  BinaryReader::readFile((heads_in[i] + ".distribs").c_str(), &distribs_loc);      
          
          ForceAssertEq(n_libs, distribs_loc.size());
          
	  for (size_t j = 0; j < n_libs; j++) {
	    const String libName = inLibNames[i][j];
	    const size_t i_lib = pairs.libraryID(libName);
            cout << Tag() << "distribution[" << setw(2) << i_lib << "]("
                 << setw(20) << libName << ") set= " << i << "," << setw(2) << j
                 << " from file" << endl;
	    distribs.at(i_lib) = distribs_loc[j];
	  }
	}
        else { // can't find file. assume gaussian.

	  for (size_t j = 0; j < n_libs; j++) {
	    const String libName = inLibNames[i][j];
	    const size_t i_lib = pairs.libraryID(libName);
            const int inv_sz_mean = pairs.getLibrarySep(i_lib);
            const int inv_sz_sd   = pairs.getLibrarySD(i_lib);
            cout << Tag() << "distribution[" << setw(2) << i_lib << "]("
                 << setw(20) << libName << ") set= " << i << "," << setw(2) << j
                 << " no file (assumed gaussian)" << endl;
            distribs.at(i_lib) = IntDistribution::gaussian(inv_sz_mean, inv_sz_sd, 4.0);
          }
        }

      }
      cout << Tag() << "merged " << distribs.size() << " distributions." << endl;
      BinaryWriter::writeFile((head_out + ".distribs").c_str(), distribs);
    }
  }






  // ---- Merge qltout files.
  
  // If the input qltout files use different reference sequences, this output
  // will be nonsensical.
  if (found_qltout) {
    cout << Tag() << "Merging qltout files" << endl;
    MergeQltout(head_out, heads_in, sizes);
  }



  
  
  // ---- Merge/Recompute kmer paths

  // If the input paths come from the same kmer numbering system, then we can
  // merge the paths simply by concatenating the paths files.  Otherwise we
  // must re-path via the slow ReadsToPathsCoreY.  The user is responsible for
  // knowing whether the kmer numbering systems are the same, and setting the
  // REPATH flag accordingly.
  

  if (!REPATH && found_paths) {
    cout << Tag() << "Merging path files" << endl;
    MergeFeudal(head_out, heads_in, paths_suffix, sizes);
  }
  
  if (REPATH) {

    cout << Tag() << "Re-creating paths files (REPATH=True)" << endl;
    vecbasevector fastb;
    for (size_t i = 0; i < n_sets; ++i)
      fastb.ReadAll((heads_in[i] + ".fastb").c_str(), true);
    
    vecKmerPath paths;
    ReadsToPathsCoreY(fastb, K, paths, head_out + ".fastb", NUM_THREADS);
    ForceAssertEq(paths.size(), size_sum);
    paths.WriteAll(head_out + paths_suffix);
  }

  // Create read tracker
 
  if (TRACK_READS) {
    cout << Tag() << "Generating read tracker" << endl;

    ReadTracker rt;

    // Add sources
    vec<uint32_t> sources(n_sets);
    for (size_t i = 0; i < n_sets; ++i)
      sources[i] = rt.AddSource(heads_in[i]);

    // Add reads
    for (size_t i = 0; i < n_sets; ++i) {
      uint32_t source_id = sources[i];
      size_t n_reads = MastervecFileObjectCount(heads_in[i] + ".fastb");
      
      for (size_t in_read_ID = 0; in_read_ID != n_reads; in_read_ID++)
	rt.AddRead(source_id, in_read_ID);
    }

    ForceAssertEq(rt.size(), size_sum);
    rt.Dump(head_out);
  }
  
  
  cout << Tag() << "Done with MergeReadSets!" << endl;
  return 0;

}


