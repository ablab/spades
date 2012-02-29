/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

const char *DOC =
"Find all useful global alignments of short queries (SEQS) against a "
"lookup table (L), fairly quickly.  Writes two files: "
"(1) OUT_PREFIX.OUT_SUFFIX contains the alignment files, in LookAlign parseable "
"format.  (2) OUT_PREFIX.minErrors.txt contains one line per read, with the number of "
"errors in the best alignment for each read, or -1 if no alignments "
"were found.  This allows you to tell if an alignment was found but not "
"written because of the MAX_ERRS or UNIQUE_ONLY filters.";

/** 
\file ShortQueryLookup.cc Find all useful global alignments of short
queries (SEQS) against a lookup table (L), fairly quickly.
*/

#include <sstream>

#include "Basevector.h"
#include "system/ParsedArgs.h"
#include "FastaFileset.h"
#include "MainTools.h"
#include "lookup/HitReceiver.h"
#include "lookup/AlignCollector.h"
#include "lookup/ImperfectLookup.h"

// Run an alignment processing chain: seqs -> look -> receiver -> aligns
// That is, look turns seqs into hits, which receiver accepts
// and turns into alignments which are passed to aligns.  The
// AlignCollectorType implements policies that determine what
// alignments are kept.  
template<class AlignCollectorType>
void LookupQueries(AlignCollectorType &aligns, lookup_table &look, 
		   const vecbasevector &seqs,
		   unsigned int MAX_INDEL_LEN, unsigned int MAX_FREQ,
		   int nPasses, int AT_START, bool DiscardReadsAtMaxFreq)
{
  if ( MAX_INDEL_LEN > 0 ) {
    GlobalHitReceiver<AlignCollectorType> 
      receiver(look, seqs, &aligns, MAX_INDEL_LEN, AT_START);
    look.FindHits(seqs, receiver, MAX_FREQ, nPasses, 0, -1, DiscardReadsAtMaxFreq);
  } else {
    ForceAssert( nPasses == 1  ||  nPasses == 2 );
    ImperfectLookup( look, seqs, aligns, nPasses == 1 ? FW : FW_OR_RC,
		     (vec<TaskTimer *>*)NULL,
		     (const vecqualvector *)NULL,
		     false /* use rc, not reverse */ );
  }
}

// Write the aligns to alignFile and the min errors to errCountFile.
// Optionally write the unique aligns to uniqueAlignFile.  The files
// are opened for appending because the command line was printed at
// the top of them earlier.
template<class AlignCollectorType>
void OutputResults(AlignCollectorType &aligns,
		   const String &errCountFile, const String &alignFile,
		   const String &uniqueAlignFile,
		   bool UNIQUE, bool READABLE, int firstRead)
{
  aligns.AdjustQueryIds(firstRead);
  {
    ofstream errStream(errCountFile.c_str(), ios_base::app);
    for (unsigned int i=0; i != aligns.size(); ++i)
      errStream << aligns.MinErrors(i) << "\n";
  }
  {
    ofstream sqltout(alignFile.c_str(), ios_base::app);
    aligns.Print(sqltout,READABLE);
  }
  if (UNIQUE) {// Print out unique results. Mostly for comparison/debugging.
    ofstream out(uniqueAlignFile.c_str(), ios_base::app); 
    for (unsigned int i=0; i != aligns.size(); ++i) {
      if (aligns.UniquelyAligned(i)) aligns.Print(out, i, READABLE);
    }
  }
}


int main( int argc, char *argv[] )
{  
  RunTime();

  BeginCommandArguments;
  CommandDoc(DOC);
  // Input parameters:
  CommandArgument_String_Doc(SEQS, "The fasta or fastb file of queries.");
  CommandArgument_String_Abbr_Doc(LOOKUP_TABLE, L,
       "The file name for a lookup table. The table's K size "
       "is used for the lookup.");
  CommandArgument_String_Abbr_OrDefault_Doc(LOOKUP_TABLE_ALT, L_ALT, "",
       "Alternate lookup table: use if LOOKUP_TABLE does not exist.");
  CommandArgument_UnsignedInt_OrDefault_Doc(K, 0, 
       "Ignored: deduced from lookup table.");
  CommandArgument_Int_OrDefault_Doc(CHUNK,50000, 
       "Process this many reads at a time. "
       "This allows us to limit the memory demands, which is "
       "critical when running the pipeline on the blades, at very "
       "little cost (for 1-chunk lookup tables.)");
  CommandArgument_Int_OrDefault_Doc(START,0, "First read to process.");
  CommandArgument_Int_OrDefault_Doc(END,-1, 
       "One past last read to align (if -1, end of file).");

  // Output parameters:
  CommandArgument_Bool_OrDefault(PRINT_UNALIGNED_READS, False);
  CommandArgument_String_Abbr_OrDefault_Doc(OUT_PREFIX, O, "",
       "Head for output file names.  If empty, use the "
       "part of SEQS before .fa.");
  CommandArgument_String_OrDefault_Doc(OUT_SUFFIX, ".sql.qltout", 
       "Extension for alignment file.");
  CommandArgument_Bool_OrDefault_Doc(READABLE, True, 
       "Write a human-readable summary of each "
       "alignment in the alignment file.");
  CommandArgument_Bool_OrDefault_Doc(UNIQUE,False, 
       "If true, output a file with all the unique "
       "alignments. This is meant mainly for comparison with "
       "UniqueUngappedLookup.");

  // Alignment search parameters:
  CommandArgument_Bool_OrDefault_Doc(FWRC, True,
       "Whether to try aligning queries in both orientations.");
  CommandArgument_UnsignedInt_OrDefault_Doc(MAX_INDEL_LEN,2, 
       "Don't look for alignments with "
       "indels larger than this.  Used to limit bandwidth "
       "for banded Smith-Waterman aligner.");
  CommandArgument_Int_OrDefault_Doc(TRIM_LEFT,0,"Align left-trimmed reads. Alignment coordinates will be returned "
				    "in terms of full read (i.e. if we trim with TRIM_LEFT=5 and trimmed "
                                    "read aligns at position Z on the ref, alignment will be saved as having "
				    "target start position Z and query start position 5, not 0). "
			    "NOTE: alignments with up to TRIM_LEFT bases overhang beyond the contig end may be returned");
  CommandArgument_Int_OrDefault_Doc(TRIM_RIGHT,0,"Trim bases from the end of each read. Alignment coordinates will be "
				    "returned in terms of full read. Alignments with up to TRIM_RIGHT bases overhang "
				    "beyond the contig end may be returned");
  CommandArgument_Int_OrDefault_Doc(MAX_FREQ,0, 
       "if positive, only use kmers occurring at most "
       "this many times in the target genome.  HACK: if maxFreq is "
       "negative, interpret it as -B, where B is the blocksize.  "
       "For each consecutive blocks of kmers, starting at the "
       "beginning of the read, take the lowest (nonzero) frequency "
       "kmer in that chunk, and put only those kmers into the query "
       "set.  That is, take one kmer from each block "
       "[0,B), [B,2B), ... .  This is a much less sensible but "
       "sometimes much faster way to search for placements in "
       "highly repetitive genomes. ");

  CommandArgument_Bool_OrDefault_Doc(IGNORE_AT_MAXFREQ,False,
       "If true, MAX_FREQ > 0 , and a single Kmer in a read has frequency "
       "higher than MAX_FREQ, then the whole read is discarded without attempting to align it. "
       "Currently works witn MAX_INDEL_LEN > 0 only");

  // Alignment filtering parameters: 
  CommandArgument_Int_OrDefault_Doc(ERR_DIFF,2, 
       "Save only alignments that have ERR_DIFF or fewer errors more than the "
       "best alignment. If negative, save all alignments. See also MAX_ERRS.");
  CommandArgument_Int_OrDefault_Doc(MAX_ERRS,4, 
       "Save only alignments that have no more than this number of errors. "
       "(plus alignments that are within ERR_DIFF of the best alignment "
       "- use ERR_DIFF=0 to exclude these extra alignments)");
  CommandArgument_Int_OrDefault_Doc(MAX_ERR_PERCENT,-1, 
       "Save only alignments that have no more "
       "than max(MAX_ERRS,MAX_ERR_PERCENT * mean_read_length) "
       "errors.");
  CommandArgument_Bool_OrDefault_Doc(UNIQUE_ONLY,False, 
       "If true, don't save or output any "
       "alignments for ambiguous reads.  This saves a great deal of "
       "time and disk space when aligning short reads to large "
       "genomes.");
  CommandArgument_Int_OrDefault_Doc(AT_START, -1, 
       "If non-negative, keep only alignments that "
       "start within AT_START bases of the beginning of a contig.");

  CommandArgument_String_OrDefault_Doc(AMB_READS_FILE,"", 
       "A list of reads which have been declared ambiguous and "
       "will not be aligned (see lookup/IdentifyAmbiguousReads.cc for details)." );
 
  EndCommandArguments;

  if ( LOOKUP_TABLE_ALT != "" && !IsRegularFile(LOOKUP_TABLE) ) 
  {    LOOKUP_TABLE = LOOKUP_TABLE_ALT;
       cout << "Using alternate lookup table.\n";    }


  if (OUT_PREFIX.empty()) {
    if ( SEQS.Contains( ".fa" ) ) OUT_PREFIX = SEQS.RevBefore( ".fa" );
    else OUT_PREFIX = SEQS;
  }

  if ( TRIM_LEFT < 0 || TRIM_RIGHT < 0 ) {
    cout << "Negative values for TRIM_LEFT or TRIM_RIGHT are not allowed" << endl;
    exit(1);
  }
  // Read the queries into seqs
  vecbasevector seqs;
  bool load_reads_by_chunks = false;
  ulonglong nreads=0; // total number of reads

  if (SEQS.Contains(".fastb")) {
    if ( MAX_ERR_PERCENT == -1 ) {
      //    seqs.ReadAll(SEQS);
      load_reads_by_chunks = true; // we can read by chunks directly from fastb, 
                                  // no need to slurp in everything and waste memory
      nreads = MastervecFileObjectCount(SEQS); // just get the total number of reads here
    } else {
      // ok, we do need all reads if we use MAX_ERR_PERCENT.
      // at least for now, until a better logic is implemented...
      seqs.ReadAll(SEQS);
      nreads = seqs.size(); 
    }
  } else {
    FastFetchReads(seqs, 0, SEQS); // can't read fasta files by chunks, have to read everything
    nreads = seqs.size(); // we loaded all the reads, so we know how many we have...
  } 

  if (-1 == END || (unsigned int)END > nreads ) END=nreads;

  PRINT3( nreads, START, END );

  const int nPasses = (FWRC ? 2 : 1);


  if ( MAX_ERR_PERCENT != -1 ) {
    // Determine mean read length and update MAX_ERRS accordingly.
    size_t sumlength=0; // total length (in bases) of all query sequences

    for (size_t i=0; i != seqs.size(); ++i) sumlength += seqs[i].size();
    // average length of a query sequence:
    double meanlength = double(sumlength)/seqs.size(); 
    MAX_ERRS= max(MAX_ERRS, int(MAX_ERR_PERCENT * meanlength / 100));
  }

  lookup_table look(LOOKUP_TABLE);
  cout << Date() << ": input read.\n";

  // Define filenames and clean up old results files if any.
  const String errCountFile = OUT_PREFIX + ".minAlignErrors.txt";
  const String alignFile = OUT_PREFIX + OUT_SUFFIX;
  const String uniqueAlignFile = OUT_PREFIX + OUT_SUFFIX + ".unique";
  Remove(errCountFile);
  Remove(alignFile);
  if (UNIQUE && !UNIQUE_ONLY)
    Remove(uniqueAlignFile);

  { // Output command header to alignFile and optionally uniqueAlignFile.
    Ofstream (alignStream, alignFile);
    if ( ! command.TheCommand( ).Contains( " NO_HEADER=True" ) && ! command.TheCommand( ).Contains( " NH=True" ) ) PrintCommandPretty(alignStream);
    if (UNIQUE && !UNIQUE_ONLY) {
      Ofstream (uniqueStream, uniqueAlignFile);
      if ( ! command.TheCommand( ).Contains( " NO_HEADER=True" ) && ! command.TheCommand( ).Contains( " NH=True" ) ) PrintCommandPretty(uniqueStream);
    }
  }
      
  vec<int> ambReads;
  if ( !AMB_READS_FILE.empty() ) {
    BinaryRead3(AMB_READS_FILE,ambReads);
    sort(ambReads.begin(),ambReads.end());
    look.SetAmbiguousReads(ambReads);
  }
  

  // Now loop over reads, from START to END-1 in groups of size CHUNK.
  // The current interval is [firstRead,lastRead).

  longlong aligned = 0;
  longlong uniquelyaligned = 0;
  for (int firstRead=START; firstRead < END; firstRead += CHUNK) {
    int lastRead = min(firstRead + CHUNK, END);
    vecbasevector smallSeqs; // Copy of relevant part of seqs
    if ( load_reads_by_chunks ) {
      // read from disk only the reads from current chunk
      smallSeqs.ReadRange(SEQS,firstRead,lastRead);
    } else {
      // copy current chunk of reads from set of all reads we pre-loaded earlier...
      smallSeqs.Append(seqs, firstRead, lastRead);
    }
    if ( TRIM_LEFT > 0 || TRIM_RIGHT > 0 ) {
      // not very efficient: causes re-writing the vectors that were
      // already copied once from the original allseds and allquals vectors
      for ( size_t i = 0 ; i < smallSeqs.size() ; i++ ) {
	smallSeqs[i].SetToSubOf(smallSeqs[i],TRIM_LEFT,smallSeqs[i].isize() - TRIM_LEFT - TRIM_RIGHT);
      }
    }
      
    cout << Date() << ", processing reads : " 
	 << firstRead << "-" << lastRead << endl;
    double startTime = WallClockTime();
    if (UNIQUE_ONLY) {
      UniqueByErrDiffAlignCollector aligns(ERR_DIFF, MAX_ERRS);
      look.SetStartRead(firstRead);
      LookupQueries(aligns, look, smallSeqs, MAX_INDEL_LEN, MAX_FREQ, 
		    nPasses, AT_START, IGNORE_AT_MAXFREQ);
      if ( TRIM_LEFT > 0 || TRIM_RIGHT > 0 ) {
	// adjust alignment start positions:
	for ( unsigned int a = 0 ; a != aligns.size(); ++a ) {
	  if ( ! aligns.UniquelyAligned(a) ) continue; // not uniquely aligned, don't bother
	  if ( aligns.Align(a).IsQueryRC() ) aligns.MutableAlign(a).a.AddToStartOnQuery(TRIM_RIGHT);
	  else aligns.MutableAlign(a).a.AddToStartOnQuery(TRIM_LEFT);
	  aligns.MutableAlign(a).query_length +=(TRIM_LEFT+TRIM_RIGHT);
	}
      }
      OutputResults(aligns, errCountFile, alignFile, uniqueAlignFile, 
		    false, READABLE, firstRead);
      aligned += AlignedCount(aligns);
      uniquelyaligned += UniquelyAlignedCount(aligns);
      if (PRINT_UNALIGNED_READS)
      {    for ( size_t j = 0; j < smallSeqs.size( ); j++ )
           {    if ( aligns.Size(j) == 0 )
                     smallSeqs[j].Print( cout, j + firstRead );    }    }
    } else {
      MaxErrDiffAlignCollector aligns(ERR_DIFF, MAX_ERRS);
      cout << " created MaxErrDiffAlignCollector" << endl;
      PRINT2( MAX_ERRS, ERR_DIFF );
      look.SetStartRead(firstRead);
      LookupQueries(aligns, look, smallSeqs, MAX_INDEL_LEN, MAX_FREQ, nPasses, 
		    AT_START, IGNORE_AT_MAXFREQ);
      aligns.Consolidate(); // This collector needs Consolidate()
      if ( TRIM_LEFT > 0 || TRIM_RIGHT > 0 ) {
	// adjust alignment start positions:
	for ( unsigned int a = 0 ; a != aligns.size(); ++a ) {
	  for ( MaxErrDiffAlignCollector::Iterator it = aligns.Begin(a) ; it != aligns.End(a) ; it++ ) {
	    if ( it->IsQueryRC() ) (*it).a.AddToStartOnQuery(TRIM_RIGHT);
	    else (*it).a.AddToStartOnQuery(TRIM_LEFT);
	  }
	}
      }
      OutputResults(aligns, errCountFile, alignFile, uniqueAlignFile, 
		    UNIQUE, READABLE, firstRead);
      PRINT2( aligned, uniquelyaligned );
      aligned += AlignedCount(aligns);
      uniquelyaligned += UniquelyAlignedCount(aligns);
      if (PRINT_UNALIGNED_READS)
      {    for ( size_t j = 0; j < smallSeqs.size( ); j++ )
           {    if ( aligns.Size(j) == 0 )
                     smallSeqs[j].Print( cout, j + firstRead );    }    }
    }
    cout << "Chunk took " << TimeSince(startTime) << endl;
  }

  // Print done message.

  const longlong N = END-START;
  cout << "Processed " << N << " reads." << endl;
  cout << "Aligned " << VALUE_AND_RATIO(3, aligned, N) << endl;
  cout << "Uniquely aligned " << VALUE_AND_RATIO(3, uniquelyaligned, N) << endl;
  String end_date = Date( );
  cout << end_date << ": done.\n";
  OfstreamMode(alignStream, alignFile, ios::app);
  alignStream << end_date << ": done.\n";
  if ( UNIQUE && !UNIQUE_ONLY ) {    
    OfstreamMode(uniqueStream, uniqueAlignFile, ios::app);
    uniqueStream << end_date << ": done.\n";    
  }

}
