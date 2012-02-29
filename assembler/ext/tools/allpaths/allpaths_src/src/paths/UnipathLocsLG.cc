///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// UnipathLocs.  Generate read locations on normal unipaths, as a 
// vec<ReadLocationLG>, where "normal" is defined by MAX_COPY_NUMBER and
// MIN_KMERS.  As compared to using all the unipaths, this saves space.
// Sort by contig, and for fixed contig, by read.  Provide a separate index by 
// read.  Files created:
// reads.unilocs.K[.rindex].

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "Intvector.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "ReadLocationLG.h"
#include "VecTemplate.h"
#include "paths/KmerPath.h"
#include "paths/PdfEntry.h"


static inline 
String Tag(String S = "ULLG") { return Date() + " (" + S + "): "; } 


int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_Int(K);
  CommandArgument_Int_OrDefault(MAX_COPY_NUMBER, 10);
  CommandArgument_Int_OrDefault(MIN_KMERS, 1);
  CommandArgument_String_OrDefault(READS, "reads");
  CommandArgument_String_OrDefault(PATHS, "paths");
  CommandArgument_String_OrDefault(UNIPATHS, "unipaths");
  CommandArgument_Bool_OrDefault_Doc(DUMP_COVERAGES, False, "Writes a report to RUN/UnipathLocs.coverages.log");
  CommandArgument_Bool_OrDefault(WRITE, True);  
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  EndCommandArguments;

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  cout << Tag() << "Using " << NUM_THREADS << " threads." << endl;

  // Set up directories.
  const String run_dir = PRE + "/" + DATA + "/" + RUN;

  vec<ReadLocationLG> ulocs;
  vec<int> ulen;
  vec<int> predicted_copyno;
  vec<bool> normal;
  size_t n_reads;
  size_t n_uni;
  {
  
    // Read in read paths and unipaths.

    const vecKmerPath paths(run_dir + "/" + READS + "." + PATHS + ".k" + ToString(K));
    const vecKmerPath paths_rc(run_dir + "/" + READS + "." + PATHS + "_rc.k" + ToString(K));
    const vecKmerPath unipaths(run_dir + "/" + READS + "." + UNIPATHS + ".k" + ToString(K));
    n_reads = paths.size();
    ForceAssertEq(n_reads, paths_rc.size());
    n_uni = unipaths.size();
    BREAD2(run_dir + "/" + READS + "." + UNIPATHS + "db.k" + ToString(K), 
           vec<tagged_rpint>, unipathsdb);
    
    cout << Tag() << "Reads to place  : " << setw(14) << n_reads << endl;
    cout << Tag() << "Unipath count   : " << setw(14) << n_uni << endl;
    
    // Compute length of each unipath.
    for (size_t i = 0; i < n_uni; i++)
      ulen.push_back(unipaths[i].KmerCount());
    
    
    {
      // Read in output of UnipathCoverage.

      VecPdfEntryVec cp;
      cp.ReadAll((run_dir + "/" + READS + "." + UNIPATHS + ".predicted_count.k" + ToString(K)).c_str());
      ForceAssertEq(cp.size(), n_uni);
      
      cout << Tag() << "Getting copy numbers." << endl;
      normal.resize(n_uni, false);
      predicted_copyno.resize(n_uni, -1);
      for (size_t i = 0; i < n_uni; i++) {
        GetMostLikelyValue(predicted_copyno[i], cp[i]);
        if (predicted_copyno[i] <= MAX_COPY_NUMBER && 
            ulen[i] >= MIN_KMERS)
          normal[i] = true;
      }
    }
    
    // Find read placements on unipaths.
    
    const size_t n_batches = NUM_THREADS;
    cout << Tag() << "Computing read locations."<< endl;
    vec< vec<ReadLocationLG> > ulocs_b(n_batches);
    #pragma omp parallel for
    for (size_t bi = 0; bi < n_batches; bi++) {
      const size_t start =  (bi    * n_reads) / n_batches;
      const size_t stop  = ((bi+1) * n_reads) / n_batches;
      
      for (size_t id = start; id < stop; id++) {
        
        for (int pass = 1; pass <= 2; pass++) {
          
          const KmerPath & p = (pass == 1 ? paths[id] : paths_rc[id]);
          vec< pair<int,int> > uo;
          
          for (int j = 0; j < p.NSegments(); j++) {
            const KmerPathInterval & I = p.Segment(j);
            vec<longlong> locs;
            Contains(unipathsdb, I, locs);
            for (int u = 0; u < locs.isize(); u++) {
              const tagged_rpint & t = unipathsdb[locs[u]];
              const int          uid = t.PathId();
              if (normal[uid]) {
                longlong offset = t.Start() - I.Start();
                for (int r = 0; r < j; r++)
                  offset += p.Segment(r).Length();
                for (int r = 0; r < t.PathPos(); r++)
                  offset -= unipaths[uid].Segment(r).Length();
                uo.push_back(make_pair(uid, int(offset)));
              }
            }
          }
          
          UniqueSort(uo);
          
          const bool orientation = (pass == 1 ? ORIENT_FW : ORIENT_RC);
          for (int t = 0; t < uo.isize(); t++) {
            const int    uid = uo[t].first;
            const int offset = uo[t].second;
            ulocs_b[bi].push_back(ReadLocationLG(id, uid, -offset, orientation)); 
          }
        }
      }
      
    }

    // Bringing everything together
 
    size_t ulocs_size = 0;
    for (size_t bi = 0; bi < n_batches; bi++)
      ulocs_size += ulocs_b[bi].size();

    ulocs.reserve(ulocs_size);
    for (size_t bi = 0; bi < n_batches; bi++)
      ulocs.append(ulocs_b[bi]);
 
  }

  // Sorting

  cout << Tag() << "Parallel sorting." << endl;
  ParallelSort(ulocs, contig_read_start_lt);
  
 
  // Generate the ulocs index.  This object allows a quick lookup, given a read_ID,
  // of all the UnipathLocs of that read.

  cout << Tag() << "Creating index." << endl;
  VecULongVec index(n_reads);
  for (size_t iii = 0; iii < ulocs.size(); ++iii)
    index[ulocs[iii].ReadId()].push_back(iii, 2, 4);
  
  // Write ulocs and index.

  if (WRITE) {
    const String filehead = run_dir + "/" + READS + ".unilocs." + ToString(K);
    cout << Tag() << "Writing locs to '" << filehead << "'." << endl;
    BinaryWrite2(filehead, ulocs);
    index.WriteAll(filehead + ".indexr");
  }
  
  // Compute some basic stats
  int unique = 0;
  int multiple = 0;
  int unplaced = 0;
  for (size_t i = 0; i < n_reads; i++) {
    const size_t placements = index[i].size();
    if      (placements == 0)  unplaced++;
    else if (placements <= 2)  unique++;
    else                       multiple++;
  }
  
  cout << Tag() << "Normal unipaths       : " << setw(14) << normal.CountValue(true)  << endl;
  cout << Tag() << "Uniquely placed reads : " << setw(14) << unique << endl;
  cout << Tag() << "Multiply placed reads : " << setw(14) << multiple << endl;
  cout << Tag() << "Unplaced reads        : " << setw(14) << unplaced << endl;
  
  // Report on the coverage of each unipath by read locs.
  if (DUMP_COVERAGES) {
    
    // Get the genome size.
    const String genome_size_file = PRE + "/" + DATA + "/genome.size";
    if (!IsRegularFile(genome_size_file)) {
      cout << Tag() 
           << "Cannot calculate coverages: file '" << genome_size_file
	   << "' does not exist.  Ignoring DUMP_COVERAGES." << endl;
    }
    else {
      const longlong genome_size = FirstLineOfFile(genome_size_file).Int();
      const double expected_coverage = double(n_reads) / genome_size;
      
      // Open the output file.
      String coverages_file = run_dir + "/UnipathLocs.coverages.out";
      cout << Tag() << "DUMP_COVERAGES: writing to '" << coverages_file << "'."<< endl;
      ofstream out(coverages_file.c_str());
      out << "Genome size = " << genome_size << endl;
      
      // Find the number of reads with a read_location on each unipath.
      vec<int> locs_per_unipath(n_uni, 0);
      for (int i = 0; i < ulocs.isize(); i++)
        locs_per_unipath[ ulocs[i].Contig() ]++;
            
      for (size_t i = 0; i < n_uni; i++) {
        int CN = predicted_copyno[i];
        
        // Find coverage of unipath by reads.
        double coverage = double(locs_per_unipath[i]) / ulen[i];
        
        out << "Unipath " << i << ":\tLength " << ulen[i] << "\tCopy # "
            << CN << "\tN reads " << locs_per_unipath[i]
            << "\t-> Coverage by read_locations: " << coverage
            << "\t(expected: " << CN << ")" << endl;
      }
      out.close();
    }
    

  }
  
  cout << Tag() << "Done!" << endl;
  return 0;
}
