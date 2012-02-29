///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// AssemblyCoverage.  Given an assembly consisting of contigs, and a
// reference sequence for the genome, assess the completeness of the
// contigs as follows.

// Find all 100 base sequences present in the reference.  For each
// such sequence, if it is present in the assembly, record it is
// covered.

// This computation takes into account multiplicity: if a 100 base
// sequence is present n times in the reference, and k times in the
// assembly, k < n, then k of the n instances are reported as covered.
// If k > n then only the n instances are reported as covered.

// WARNING: scaffold breaks are determined by nnnnnnnnnnnnnnnnnnn sequences.

// Now accepts efasta as input.

// Note that the creation of kmers involving ambiguities may introduce the
// appearance of overcoverage.  Conversely, note that kmer multiplicity is capped.

#include "Basevector.h"
#include "Bitvector.h"
#include "Fastavector.h" // fastavector
#include "FastIfstream.h" // fast_ifstream
#include "FetchReads.h"
#include "FetchReadsAmb.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "SeqInterval.h"
#include "efasta/EfastaTools.h"
#include "kmers/KmerParcelsTools.h" // CoupledParcelKmers
#include "math/Functions.h"
#include "reporting/PerfStat.h"
#include "util/RunCommand.h"

int main(int argc, char **argv)
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String_Doc(ASSEMBLY, 
    "fasta file for assembly; Ns are treated as gaps");
  CommandArgument_String_Doc(REF, "fasta file for reference");
  CommandArgument_String_Doc(HEAD,
    "directory for temporary files (eg, parcel data)");
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_CONTIG, 1000,
    "imposed on both ASSEMBLY and REF");
  CommandArgument_Int_OrDefault(K, 100);
  CommandArgument_Int_OrDefault(MAX_KMER_MULTIPLICITY, 10);
  CommandArgument_Int_OrDefault(MAX_PER, 100);
  CommandArgument_Bool_OrDefault_Doc(ARCHIVE, False,
    "send output to <HEAD>.log rather than cout");
  CommandArgument_Bool_OrDefault_Doc(OVERWRITE, True,
    "if false, it does not regenerate data if it already exists");
  CommandArgument_Bool_OrDefault_Doc(CLEAN_UP, True,
    "if true, remove intermediate files (kmer parcels) at the end");
  EndCommandArguments;
  
  // Dir and file names.
  String parcels_ref_dir = HEAD + "/reference";
  String parcels_ass_dir = HEAD + "/assembly";

  String splitref_fastb_file = HEAD + "/splitref.fastb";
  String splitass_fastb_file = HEAD + "/splitass.fastb";
  String splitref_map_file = HEAD + "/splitref.map";
  String kwin_coverage_file = HEAD + "/coverage.kwin";

  // Thread control

  NUM_THREADS = configNumThreads(NUM_THREADS);

  Mkpath(HEAD);

  String log_file = HEAD + "/AssemblyCoverage.log";
  ostream *pout = 0;
  if (ARCHIVE) {
    pout = (ostream*) new ofstream(log_file.c_str());
    cout << "\nSending output to " << log_file << "\n" << endl;
    PrintCommandPretty(*pout);
  }
  ostream &out = *(ARCHIVE ? pout : (ostream*) &cout);
  
  // Load the reference (split at ambiguous bases).
  out << Date() << ": loading reference genome" << endl;
  vecbasevector ref;
  vecbitvector is_covered;   // WARNING: in kmer space!
  vec<seq_interval> cgmap;
  vec<int> rklens;   // WARNING: in kmer space!
  {
    vecbasevector ref0;
    vecbitvector ref_amb;
    FetchReads(ref0, 0, REF);
    FetchReadsAmb(ref_amb, REF, AMB_EQ_ANY_AMBIGUOUS);
    rklens.resize(ref0.size(), 0);

    Mimic(ref0, is_covered, False);
    for (int ii=0; ii<(int)is_covered.size(); ii++)
      is_covered[ii].resize(Max(0, (int)is_covered[ii].size() - (K-1)));

    out << Date() << ": Breaking reference genome on ambiguities" << endl;
    for (size_t i = 0; i < ref0.size(); i++) {
      rklens[i] = ref0[i].size() - (K-1);
      vec<char> C;
      for (int j = 0; j < ref0[i].isize(); j++) {
	if (ref_amb[i][j] && C.nonempty()) {
	  basevector b(C.size()); 
	  for (int k = 0; k < C.isize(); k++) b.Set(k, C[k]);
	  if (b.size() >= MIN_CONTIG) {
	    int interval_id = ref.size();
	    int begin = j - C.isize();
	    int end = j;
	    cgmap.push_back(seq_interval(interval_id, i, begin, end));
	    ref.push_back_reserve(b);
	  }
	  C.clear();
	}
	if (!ref_amb[i][j]) C.push_back(ref0[i][j]);
      }
      if (C.nonempty()) {
	basevector b(C.size());
	for (int k = 0; k < C.isize(); k++) b.Set(k, C[k]);
	if (b.size() >= MIN_CONTIG) {
	  int interval_id = ref.size();
	  int begin = ref0[i].isize() - C.isize();
	  int end = ref0[i].isize();
	  cgmap.push_back(seq_interval(interval_id, i, begin, end));
	  ref.push_back_reserve(b);
	}
      }
    }
  }
  
  // Save split reference fastb and map.
  out << Date() << ": saving split reference to '" << splitref_fastb_file << "'." << endl;
  ref.WriteAll(splitref_fastb_file);
  
  ofstream siout(splitref_map_file.c_str());
  for (int ii=0; ii<cgmap.isize(); ii++)
    siout << cgmap[ii] << "\n";
  siout.close();
  
     // Load assembly (split at ambiguous bases).
     
     vecbasevector assembly;
     if ( OVERWRITE || ( !IsRegularFile(splitass_fastb_file) ) ) 
       {    out << Date( ) << ": loading assembly from '" << ASSEMBLY << "'."<< endl;
          Bool EFASTA = ASSEMBLY.Contains( ".efasta", -1 );
          if ( !EFASTA )
          {    vec<fastavector> afv;
               LoadFromFastaFile(ASSEMBLY, afv);
               out << Date() << ": splitting " << afv.size() 
                    << " fastavectors (.=100 fastavectors)" << endl;
               for (size_t ii=0; ii<afv.size(); ii++) 
               {    if (ii % 100 == 0) Dot(out, ii / 100);
                    vec<fastavector> chunks = afv[ii].SplitOnGaps( );
                    for (size_t jj=0; jj<chunks.size(); jj++) 
                    {    if (chunks[jj].size() < MIN_CONTIG) continue;
	                 vecbvec bases 
                              = chunks[jj].AllKmers(K, MAX_KMER_MULTIPLICITY);
	                 for (size_t kk=0; kk<bases.size(); kk++)
	                      assembly.push_back_reserve(bases[kk]);    }    }
               out << endl;

               // Commented out, as the split assembly can be huge and very 
               // slow to save.
               // out << Date() << ": saving split assembly" << endl;
               // assembly.WriteAll(splitass_fastb_file);
                    }
          else
          {    vec<efasta> contigs, scaffolds;
               LoadEfastaIntoStrings( ASSEMBLY, scaffolds );
	       cout << Date() << ": " << scaffolds.size() << " efasta scaffolds loaded." << endl;
               vec<superb> superbs;
               SplitEfastaIntoContigs( scaffolds, contigs, superbs );    
	       cout << Date() << ": " << contigs.size() << " contigs from scaffolds." << endl;
               for ( size_t i = 0; i < contigs.size( ); i++ )
               {    if ( contigs[i].size( ) < MIN_CONTIG ) continue;
                    vecbasevector bases;
                    contigs[i].GetKmers( K, bases, MAX_PER );
                    for (size_t kk=0; kk<bases.size(); kk++)
                         assembly.push_back_reserve(bases[kk]);    }    }    }
     else 
     {    out << Date() << ": loading split assembly" << endl;
          assembly.ReadAll(splitass_fastb_file);    }
  
  // Find shared kmers between the reference and assembly.

  out << Date( ) << ": loading assembly from '" << ASSEMBLY << "'."<< endl;
  vecbasevector shared_kmers;
  vec<int> ref_freq, ass_freq;
  {
    bool build_parcels = OVERWRITE;

    KmerParcelsDiskStore parcels_ass(K, parcels_ass_dir);
    KmerParcelsDiskStore parcels_ref(K, parcels_ref_dir);

    if (!parcels_ass.ExistOnDisk()) build_parcels = true;
    if (!parcels_ref.ExistOnDisk()) build_parcels = true;

    if (build_parcels) {
      out << Date() << ": Building coupled " << K << "-mer parcels." << endl;

      parcels_ref.SetKeepOnDisk(!CLEAN_UP);
      parcels_ass.SetKeepOnDisk(!CLEAN_UP);
    
      KmerParcelsCoupledBuilder builder(assembly, ref, 
                                        parcels_ass, parcels_ref, 
                                        NUM_THREADS);
      
      builder.Build();
    }
    else {
      out << Date() << ": parcels on disk, no need to build." << endl;
    }

    size_t n_parcels = parcels_ref.GetNumParcels();


    // Reserve memory for the data structures we need.
    // ribeiro: probably an overestimation, since it is based on the
    //          number of of batches in the reference, and not on the
    //          number of matching batches.
    // ribeiro: do we really need to reserve?  Is it such a
    //          performance drain not to?  Just asking.
    size_t n_batches_ref = parcels_ref.GetTotalNumKmerBatches();
    size_t n_bases = n_batches_ref * K;
    shared_kmers.Reserve(n_batches_ref + n_bases/16, n_batches_ref);
    ref_freq.reserve(n_batches_ref);
    ass_freq.reserve(n_batches_ref);
    
    // Open up the intersection KmerParcels and process the batches,
    // finding the (shared) basevector and the frequency in each dataset.
    out << Date() << ": parsing parcels" << endl;
    for (size_t parcel_ID = 0; parcel_ID < n_parcels; parcel_ID++) {
      KmerParcelReader parcel_ref(parcels_ref, parcel_ID);
      KmerParcelReader parcel_ass(parcels_ass, parcel_ID);
      
      while (GetNextMatchingKmerBatches(parcel_ref, parcel_ass)) {
        const KmerBatch & batch_ref = parcel_ref.CurrentKmerBatch();
        const KmerBatch & batch_ass = parcel_ass.CurrentKmerBatch();
        
        BaseVec kmer;
        batch_ref.GetKmerBaseVec(&kmer);
	shared_kmers.push_back(kmer);
	
	ref_freq.push_back(batch_ref.GetKmerFreq());
	ass_freq.push_back(batch_ass.GetKmerFreq());

	// Record shared kmer as covered on reference.
	for (int ii = 0; ii < Min(batch_ref.isize(), 10); ii++) {
	  int id = batch_ref[ii].GetReadID();
	  int pos = batch_ref[ii].GetUnsignedPos();
	  int super_id = cgmap[id].SeqId();
	  int super_offset = cgmap[id].Begin();
	  is_covered[super_id].Set(pos + super_offset, true);
	}
      }
    }
  }
  
  // Turn bitvector into windows of coverage (in kmer space).
  vec<seq_interval> ckwin;
  for (int ii=0; ii<(int)is_covered.size(); ii++) {
    int currsize = (int)is_covered[ii].size();
    int pos = 0;
    while (pos < currsize) {
      while (pos < currsize && ! is_covered[ii][pos]) pos++;
      int begin = pos;
      while (pos < currsize && is_covered[ii][pos]) pos++;
      int end = pos;
      if (begin < end) {
	int interval_id = ckwin.size();
	ckwin.push_back(seq_interval(interval_id, ii, begin, end));
      }
    }
  }

  // Save windows of coverage.
  out << Date() << ": saving windows of coverage" << endl;
  {
    ofstream out(kwin_coverage_file.c_str());
    for (int ii=0; ii<ckwin.isize(); ii++)
      out << ckwin[ii] << "\n";
    out.close();
  }
  
  // Compute coverage.
  out << Date() << ": Computing coverage" << endl;
  
  // Find the "shared frequency", which is the minimum of the two freqs.
  size_t n_shared_kmers = shared_kmers.size();
  vec<int> shared_freqs(n_shared_kmers);
  for (size_t i = 0; i < n_shared_kmers; i++)
    shared_freqs[i] = Min(ref_freq[i], ass_freq[i]);
  
  // Count the total number of kmers.
  longlong total = 0, GC_total = 0;
  for (size_t i = 0; i < ref.size(); i++)
    if (ref[i].isize() >= K) {
      total += ref[i].size() - K + 1;
      GC_total += GcBases(ref[i]);
    }
  
  longlong covered = BigSum(shared_freqs);
  longlong uncovered = total - covered;
  longlong overcovered = 0;
  for (size_t i = 0; i < n_shared_kmers; i++)
    overcovered += Max(0, ass_freq[i] - ref_freq[i]);
  
  // Count GC content in the covered bases.
  longlong GC_covered = 0, GC_uncovered = 0;
  for (size_t i = 0; i < n_shared_kmers; i++)
    GC_covered += GcBases(shared_kmers[i]) * shared_freqs[i];
  
  // Look for unique kmers in the reference and count their coverage.
  longlong covered_unique = 0, uncovered_unique = 0;
  for (size_t i = 0; i < n_shared_kmers; i++)
    if (ref_freq[i] == 1) covered_unique++;
  
  out << endl;
  out << "REPORT\n" << endl;
  out << setiosflags(ios::fixed) << setprecision(1) 
       << 100.0 * double(covered)/double(total) << resetiosflags(ios::fixed) 
       << "% covered" << endl;
  if ( covered == total )
  {    PerfStat::log( ) << PerfStat( "frac_covered", 
            "% of genome covered by contigs", 100 );    }
  else
  {    PerfStat::log( ) << std::fixed << std::setprecision(1)
            << PerfStat( "frac_covered", "% of genome covered by contigs",
            100.0 * double(covered) / double(total) );    }
  out << PERCENT_RATIO(3, overcovered, total) << " duplicated" << endl;
  PerfStat::log( ) << std::fixed << std::setprecision(2)
       << PerfStat( "frac_dup", "% of duplication in contigs",
       100.0 * double(overcovered) / double(total) );
  out << "GC content of covered bases: "
       << PERCENT_RATIO(3, GC_covered, covered*longlong(K)) << endl;
  out << "GC content of total bases: "
       << PERCENT_RATIO(3, GC_total, total) << endl;
  out << "unique fraction in covered = "
       << PERCENT_RATIO(3, covered_unique, covered) << endl;
  //out << "unique fraction in uncovered = "
  //     << PERCENT_RATIO(3, uncovered_unique, uncovered) << endl;
  out << endl;


  // Done.
  cout << Date() << ": done" << endl;
  if (pout) {
    out << Date() << ": done" << endl;
    delete(pout);
  }
  
}

