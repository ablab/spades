///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* CacheReadsMerge
 * 
 * A module to take a set of reads and library info and 
 * convert it into a form usable by Arachne pipelines such as RunAllPathsLG.
 *
 * WARNING: This module has not been tested with unpaired lanes.
 *
 * ADVICE: This module is hard to run by hand.  It expects a lot of info from
 *         the command line which is much easily provided by scripts.
 *
 * Command-line arguments:
 *
 *   CACHE_DIR: Root directory for the pre-processed input data.
 * 
 *   GROUP_HEADS: list of all the file group heads to merge.
 *
 *   GROUP_READ_SIZES: list of mean read sizes in the groups.
 *
 *   GROUP_LIB_NAMES: list of library names corresponding to the GROUP_HEADS.
 *
 *   GROUP_LIB_PAIRED: list of 1s or 0s for paired or unpaired libraries, respectively.
 *
 *   GROUP_LIB_SIZES: list of library sizes.
 *
 *   GROUP_LIB_STDDEVS: list of library standard deviations.
 *
 *   OUT_HEAD: prefix for the output files.
 *
 *
 * The following files are created:
 *
 *   <OUT_HEAD>.fastb
 *   <OUT_HEAD>.qualb
 *   <OUT_HEAD>.pairs
 *   <OUT_HEAD>.source.txt
 *
 * NOTE: Most of the code was pilfered from Josh's solexa/FetchIluminaReads.cc
 * 
 *
 * 2010-05-28    Filipe Ribeiro     <crdhelp@broadinstitute.org>
 *
 ******************************************************************************/
#include "MainTools.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "String.h"
#include "lookup/LibInfo.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/VirtualMasterVec.h"
#include "random/Shuffle.h"

typedef VirtualMasterVec<BaseVec> VBaseVecVec;
typedef VirtualMasterVec<QualVec> VQualVecVec;


static inline 
String Tag(String S = "CRM") { return Date() + " (" + S + "): "; } 



void vec_to_text_file(const vec<size_t> & v,
                      const String & fn)
{
  ofstream os(fn.c_str());
  const size_t n = v.size();
  os << n << endl;
  for (size_t i = 0; i != n; i++)
    os << v[i] << endl;
  os.close();
}





vec<size_t> get_shuffled_read_indices(const size_t n_in_reads, 
                                      const size_t n_out_reads, 
                                      const int LIB_PAIRED,
                                      const unsigned int SEED)
{
  vec<size_t> ir_shuffle;

  if (LIB_PAIRED) {
    const size_t n_pairs_in  = n_in_reads / 2;
    const size_t n_pairs_out = n_out_reads / 2;
    
    // ---- get shuffled pair indices (assumes interleaved reads) 
    vec<uint64_t> ip_shuffle;
    Shuffle64(n_pairs_in, ip_shuffle, SEED);
    ip_shuffle.resize(n_pairs_out);
    sort(ip_shuffle.begin(), ip_shuffle.end());
    
    // ---- get shuffled read indices (assumes interleaved reads) 
    for (size_t j = 0; j != n_pairs_out; j++) {
      const size_t ip = ip_shuffle[j];
      ir_shuffle.push_back(2 * ip);
      ir_shuffle.push_back(2 * ip + 1);
    }
  }
  else { // library is not paired
    // ---- get shuffled reads indices 
    Shuffle64(n_in_reads, ir_shuffle, SEED);
    ir_shuffle.resize(n_out_reads);
    sort(ir_shuffle.begin(), ir_shuffle.end());
  }

  return ir_shuffle;
}








void merge_files(const vec<String> & in_bases_fns, 
                 const vec<String> & in_quals_fns, 
                 const vec<String> & FILE_HEADS, 
                 const vec<int> & LIB_PAIRED,
                 const vec<size_t> & n_in_reads,
                 const vec<size_t> & n_out_reads,
                 const String & OUT_HEAD, 
                 const unsigned int SEED,
                 const bool do_quals) 
{
  const size_t n_files = in_bases_fns.size();
  const String out_bases_fn = OUT_HEAD + ".fastb";
  const String out_quals_fn = OUT_HEAD + ".qualb";

  bool all_reads = true;
  for (size_t i = 0; i != n_files; i++)
    all_reads = (all_reads && (n_in_reads[i] == n_out_reads[i]));


  if (all_reads) { // ---- merge all reads
    
    if (n_files == 1) {
      Cp(in_bases_fns[0], out_bases_fn); 
      if (do_quals)
        Cp(in_quals_fns[0], out_quals_fn);
    } 
    else {
      MergeMastervecs(out_bases_fn, in_bases_fns);
      if (do_quals)
        MergeMastervecs(out_quals_fn, in_quals_fns);
    }
  }
  else {  // ---- merge random fractions of reads

    // ---- declare writers
    IncrementalWriter<BaseVec> bv_out(out_bases_fn.c_str());
    IncrementalWriter<QualVec> qv_out(out_quals_fn.c_str());

    for (size_t i = 0; i != n_files; i++) {
      cout << Tag() << "From '" << FILE_HEADS[i] << ".fastb'" 
           << " [" << ((LIB_PAIRED[i]) ? "paired" : "unpaired") << "]" 
           << ": " << n_out_reads[i] << " out of " << n_in_reads[i] << " reads"
           << " (" << (100.0 * n_out_reads[i] / n_in_reads[i]) << " %)" 
           << endl;

      // ---- get shuffled read indices (paired or not)
      const vec<size_t> ir_shuffle = get_shuffled_read_indices(n_in_reads[i], n_out_reads[i], 
                                                               LIB_PAIRED[i], SEED);
      
      // ---- output indices of selected reads
      vec_to_text_file(ir_shuffle, OUT_HEAD + "." + FILE_HEADS[i] + ".i_reads");
        


      // ---- load selected read subset and merge
      const size_t nr = ir_shuffle.size();

      VBaseVecVec bv_in(in_bases_fns[i].c_str());
      for (size_t j = 0; j != nr; j++)
        bv_out.add(bv_in[ir_shuffle[j]]);

      if (do_quals) {
        VQualVecVec qv_in(in_quals_fns[i].c_str());
        for (size_t j = 0; j != nr; j++)
          qv_out.add(qv_in[ir_shuffle[j]]);
      }
    }

    bv_out.close();
    qv_out.close();
  }
  
}





int main(int argc, char **argv)
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String(CACHE_DIR);
  CommandArgument_StringSet(GROUP_HEADS);
  CommandArgument_DoubleSet(GROUP_READ_SIZES);
  CommandArgument_DoubleSet(GROUP_FRACS);
  CommandArgument_StringSet(GROUP_LIB_NAMES);
  CommandArgument_IntSet(GROUP_LIB_PAIRED);
  CommandArgument_IntSet(GROUP_LIB_SIZES);
  CommandArgument_IntSet(GROUP_LIB_STDDEVS);
  CommandArgument_UnsignedInt_OrDefault( SEED, 666 );
  CommandArgument_String(OUT_HEAD);
  EndCommandArguments;
  

 
  // Sanity checks: Determine if any of the input arguments are clearly invalid.
  if (!IsDirectory(CACHE_DIR))
    FatalErr("CACHE_DIR=" << CACHE_DIR << " does not exist.");

  if (!OUT_HEAD.Contains("/")) 
    OUT_HEAD = String(getenv("PWD")) + "/" + OUT_HEAD;

  const size_t n_groups = GROUP_HEADS.size();

  ForceAssertEq(GROUP_LIB_NAMES.size(), n_groups);
  ForceAssertEq(GROUP_LIB_PAIRED.size(), n_groups);
  ForceAssertEq(GROUP_LIB_SIZES.size(), n_groups);
  ForceAssertEq(GROUP_LIB_STDDEVS.size(), n_groups);
  ForceAssertEq(GROUP_READ_SIZES.size(), n_groups);
  ForceAssertEq(GROUP_FRACS.size(), n_groups);
  for (size_t i = 0; i < n_groups; i++)
    ForceAssert(GROUP_FRACS[i] >= 0.0 && GROUP_FRACS[i] <= 1.0);




  // ---- Collect some data and verify consistency of reads and quals

  cout << Tag() << "Starting merge of fastb/qualb files from all groups..." << endl;

  vec<String> in_bases_fns;
  vec<String> in_quals_fns;
  vec<size_t> n_in_reads;
  vec<size_t> n_out_reads;

  bool do_quals = true;

  for (size_t i = 0; i < n_groups; i++) {
    const String in_head = CACHE_DIR + "/" + GROUP_HEADS[i];
    in_bases_fns.push_back(in_head + ".fastb");
    in_quals_fns.push_back(in_head + ".qualb");

    const size_t n_r = MastervecFileObjectCount(in_bases_fns[i]);
    
    if (IsRegularFile(in_quals_fns[i]))
      ForceAssertEq(n_r, MastervecFileObjectCount(in_quals_fns[i]));
    else 
      do_quals = false;

    n_in_reads.push_back(n_r);
    n_out_reads.push_back(GROUP_LIB_PAIRED[i] ? 
                          2 * size_t((GROUP_FRACS[i] * n_r) / 2) : 
                          size_t((GROUP_FRACS[i] * n_r)));
  }





 
   
  
  // ---- Actually do the merging 

  cout << Tag() << "Merging read set into '" << OUT_HEAD << ".{fastb,qualb}'." << endl;

  merge_files(in_bases_fns, in_quals_fns, GROUP_HEADS, 
              GROUP_LIB_PAIRED, n_in_reads, n_out_reads,
              OUT_HEAD, SEED, do_quals);

  const size_t n_out_reads_total = Sum(n_out_reads);
  const size_t n_in_reads_total = Sum(n_in_reads);
  cout << Tag() << setw(12) << n_in_reads_total << "  reads in." << endl;
  cout << Tag() << setw(12) << n_out_reads_total << "  reads out." << endl;
  cout << Tag() << setw(12) << fixed << setprecision(1)
       << (100.0 * n_out_reads_total / n_in_reads_total) << "  % of reads included." 
       << endl;






  
  // ---- Create a PairsManager


  cout << Tag() << "Creating PairsManager." << endl;
  PairsManager pairs(n_out_reads_total);
  {
    map<String, float> n_sz_sum;
    map<String, float> n_sum;

    for (size_t i = 0; i < n_groups; i++) {
      const String lib_name = GROUP_LIB_NAMES[i];
      n_sz_sum[lib_name] += n_out_reads[i] * GROUP_READ_SIZES[i];
      n_sum[lib_name] += n_out_reads[i];
    }

    for (size_t i = 0; i < n_groups; i++) {
      const String lib_name = GROUP_LIB_NAMES[i];
      if (pairs.libraryID(lib_name) < 0) { // library not in pairs yet
        const int sz_mean = n_sz_sum[lib_name] / n_sum[lib_name];
        const int sep_mean = GROUP_LIB_SIZES[i] - 2 * sz_mean;
        pairs.addLibrary(sep_mean, GROUP_LIB_STDDEVS[i], GROUP_LIB_NAMES[i]);
      }
    }
  }
  



  // For each lane of paired reads, find the read pair separation/stdev,
  // and then add a set of read pairings to the PairsManager.
  // Here we implicitly assume that paired reads are interleaved in the 
  // fastb: that is, reads 0 and 1 are paired, 2 and 3 are paired, etc.
  
  for (size_t i = 0, i0_read = 0; i < n_groups; i0_read += n_out_reads[i++]) {
    
    if (GROUP_LIB_PAIRED[i]) {  // skip unpaired libraries

      const int i_lib = pairs.libraryID(GROUP_LIB_NAMES[i]);  
      ForceAssertGe(i_lib, 0); 
      
      for (size_t j = 0; j < n_out_reads[i]; j += 2)  // assume interleaving
	pairs.addPairToLib(i0_read + j, i0_read + j + 1, i_lib);
    }
  }



  // ---- Save PairsManager.
  pairs.Write(OUT_HEAD + ".pairs");




  // ---- Output some statistics
  cout << Tag() << "Library statistics:" << endl << endl;
  pairs.printLibraryStats(cout);
  cout << endl;


  // ---- create a text file indicating the source of this dataset.
  ofstream source((OUT_HEAD + ".source.txt").c_str());
  command.PrintTheCommandPretty(source);

  // ---- write pairs stats to the summary file
  pairs.printLibraryStats(source);
  source.close();
  
  
  cout << Tag() << ": Done!" << endl;
}
