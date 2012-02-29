/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Evaluate error correction using reference, not true reads.

const char *DOC =
"EvalCorrectedReads. Evaluates success of error correction and suspicious "
"read removal by perfectly aligning these reads back to a reference. "
"Does not require truth data";

#include "Basevector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "lookup/PerfectLookup.h"
#include "lookup/LookAlign.h"


String Fraction(longlong numerator, longlong denominator) {
  if (denominator == 0)
    return "N/A";
  return ToString(double(numerator)/denominator, 3);
}


// Add aligned pair to coverage map. For multiple paired aligments, add fractional
// coverage - 1/alignments to each base position in the map.
void AddCoverage(const vec<look_align>& aligns, const vec<pair<longlong, longlong> >& valid,
		 vec<vec< double> >& covered) {
  if (valid.empty())
    return;
  longlong alignments = valid.size();
  double fraction = 1.0/alignments;
  for (longlong i = 0; i < alignments; i++) {
    for (int read = 0; read < 2; read++) {  // 1st or 2nd read in pair
      const look_align& la = 
	(read == 0 ? aligns[valid[i].first] : aligns[valid[i].second]);
      int contig = la.TargetId();
      int start = la.StartOnTarget();
      int end = la.EndOnTarget();
      for (int pos = start; pos <= end; pos++)
	covered[contig][pos] += fraction;
    }   
  }
}


// Write base coverage statistics
void CoverageStatistics(vec<vec< double> >& covered, int coverage) {
  longlong low00 = 0;
  longlong low25 = 0;
  longlong low50 = 0;
  longlong low75 = 0;
  size_t total = SizeSum(covered);
  for (longlong i = 0; i < covered.isize(); i++) {
    low00 += covered[i].CountValue(0);
    low25 += count_if(covered[i].begin(), covered[i].end(),
		      bind2nd(less<double>(), coverage/4));
    low50 += count_if(covered[i].begin(), covered[i].end(),
		      bind2nd(less<double>(), coverage/2));
    low75 += count_if(covered[i].begin(), covered[i].end(),
		      bind2nd(less<double>(), coverage*3/4));
  }
  cout << "\nBases covered by paired reads at...:\n";
  cout << "  0-25% of expected coverage : " << low25 
       <<"  (" << Fraction(low25, total) << ")\n";
  cout << "  0-50% of expected coverage : " << low50
       <<"  (" << Fraction(low50, total) << ")\n";
  cout << "  0-75% of expected coverage : " << low75 
       <<"  (" << Fraction(low75, total) << ")\n";
  cout << "  75%+ of expected coverage  : " << (total - low75) 
       <<"  (" << Fraction(total - low75, total) << ")\n";
  cout << "  Uncovered bases            : " << low00 
       <<"  (" << Fraction(low00, total) << ")\n";
}

// Evaluate read pair given set of possible alignments of each read
// Return true if at least one sensible alignment is found
// Number of alignments found is returned in alignmentCount
// Multiple alignments are allowed, and each possible pair alignment
// is returned as a pair of alignment index positions in vec valid.
// A sensible alignment is one that has an actual separation within 
// max*expectedSd of the expected separation.
Bool EvaluatePair(const vec<look_align>& aligns, const vec<int>& indexReadA,
		  const vec<int>& indexReadB, vec<pair<longlong, longlong> >& valid,
		  int expectedSep, int expectedSd,
		  int& alignmentCount, double max = 3) {

  if (indexReadA.empty() || indexReadB.empty()) // need both to align
    return False;

  valid.clear();

  alignmentCount = 0;
  for (longlong iA = 0; iA < indexReadA.isize(); iA++) {
    for (longlong iB = 0; iB < indexReadB.isize(); iB++) {
      const look_align& laA = aligns[indexReadA[iA]];
      const look_align& laB = aligns[indexReadB[iB]];
      int startA = laA.StartOnTarget();
      int startB = laB.StartOnTarget();
      if (laA.rc1 == laB.rc1) // incorrect orientation
	continue;
      if (laA.TargetId() != laB.TargetId())
	continue;  // land on different reference contigs
      int actualSep;
      if (laA.rc1) 
	actualSep = startA-startB-laB.query_length;
      else 
	actualSep = startB-startA-laA.query_length;
      if (abs(actualSep - expectedSep) < max * expectedSd) {
	alignmentCount++;
	valid.push_back(make_pair(indexReadA[iA], indexReadB[iB]));
      }
    }
  }
  return (alignmentCount);
}

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_Int_OrDefault_Doc(LOOKUP_K, 12,
    "Kmer size of lookup table");
  CommandArgument_Int_OrDefault_Doc(MAX_NO_SD, 3,
    "Maximum number of standard deviation allowed from mean separation for a "
    "pair to be still consider valid");
  CommandArgument_String_OrDefault_Doc(READS_ORIG, "",
    "Uncorrected reads.");
  CommandArgument_String_OrDefault_Doc(READS_IN, "reads_nosusp",
    "Reads that we need to find locations on the genome for.");
  CommandArgument_String_OrDefault(GENOME, "genome");
//   CommandArgument_Bool_OrDefault_Doc(UNIQUE_ONLY, False,
//     "If True, only use uniquely place reads." );
  EndCommandArguments;

  // Set up directories.
  String data_dir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;

  size_t genomeSize = 0;
  size_t genomeContigCount = 0;
  vec<unsigned int> genomeContigSizes;

  // Determine genome size statistics
  { 
    vecbasevector genome(data_dir + "/" + GENOME + ".fastb");
    genomeSize = genome.sumSizes();
    genomeContigCount = genome.size();
    genomeContigSizes.resize( genomeContigCount );
    for (size_t i = 0; i < genomeContigCount; i++)
      genomeContigSizes[i] = genome[i].size( );
  }

  // Determine original number of reads and pairs
  longlong originalReads = 0;
  longlong originalPairs = 0;
  if (READS_ORIG != "") {
    originalReads = MastervecFileObjectCount( run_dir + "/" + READS_ORIG + ".fastb" );
    PairsManager orig_pairs( run_dir + "/" + READS_ORIG + ".pairs" );
    originalPairs = orig_pairs.nPairs();
  }

  // Load corrected reads
  vecbasevector reads(run_dir + "/" + READS_IN + ".fastb");
  longlong n_reads = reads.size();

  // Load pairing information
  PairsManager pairs( run_dir + "/" + READS_IN + ".pairs");
  size_t n_pairs = pairs.nPairs();

  // Compute estimated coverage by paired reads
  longlong estimatedPairCoverage = 0;
  for (size_t i = 0; i < n_pairs; i++ ) {
    longlong id1 = pairs.ID1(i);
    longlong id2 = pairs.ID2(i);
    estimatedPairCoverage += reads[id1].size() + reads[id2].size();
  }
  PRINT3( estimatedPairCoverage, genomeSize, estimatedPairCoverage / genomeSize );
  
  cout << "Aligning " << n_reads << " reads to reference" << endl;

  // Perform Alignment
  vec<look_align> lookAligns;
  PerfectLookup(LOOKUP_K, reads, data_dir + "/" + GENOME + ".lookup", lookAligns, FW_OR_RC);

  cout <<  "Found " << lookAligns.size() << " read alignments " << endl;

  // Build alignment index
  vec<vec<int> > index(n_reads);
  BuildLookAlignsIndex(lookAligns, index, n_reads);

  // Count number of aligned reads
  longlong alignedReads = 0;
  for (longlong i = 0; i < index.isize(); i++ ) 
    if (index[i].nonempty())
      alignedReads++;


  longlong uniquelyAlignedPairs = 0;
  longlong multiplyAlignedPairs = 0;
  longlong unalignedPairs = 0;
  longlong pairCoverage = 0;

  vec<pair<longlong, longlong> > valid;

  vec< vec< double> > covered(genomeContigCount);
  for (size_t i = 0; i < genomeContigCount; i++)
    covered[i] = vec<double>(genomeContigSizes[i],0);

  // Evaluate paired alignments
  for (size_t i = 0; i < n_pairs; i++ ) {
    longlong id1 = pairs.ID1(i);
    longlong id2 = pairs.ID2(i);
    int aligned;
    if (EvaluatePair(lookAligns, index[id1], index[id2], valid,
		     pairs.sep(i), pairs.sd(i), aligned, MAX_NO_SD) ){
      if (aligned == 1 )
	uniquelyAlignedPairs++;
      else 
	multiplyAlignedPairs++;
      pairCoverage += reads[id1].size() + reads[id2].size();
      AddCoverage(lookAligns, valid, covered);
    } else
      unalignedPairs++;
  }

  longlong alignedPairs = uniquelyAlignedPairs + multiplyAlignedPairs;

  // Write out some statistics

  cout << "\nGenome Size  : " << genomeSize << "\n";

  cout << "\nReads:\n";
  cout << "  Original   : " << originalReads << "\n";
  cout << "  Corrected  : " << n_reads 
       <<"  (" << Fraction(n_reads, originalReads) << " of original)\n";
  cout << "  Aligned    : " << alignedReads
       <<"  (" << Fraction(alignedReads, n_reads) << " of corrected)\n";
  cout << "  Unaligned  : " << n_reads - alignedReads
       <<"  (" << Fraction(n_reads - alignedReads, n_reads) << " of corrected)\n";
  
  cout << "\nPaired Reads:\n";
  cout << "  Original   : " << originalPairs << "\n";
  cout << "  Corrected  : " << n_pairs
       <<"  (" << Fraction(n_pairs, originalPairs) << " of original)\n";
  cout << "  Aligned    : " << alignedPairs 
       <<"  (" << Fraction(alignedPairs, n_pairs) << " of corrected)\n";
  cout << "    Uniquely : " << uniquelyAlignedPairs
       <<"  (" << Fraction(uniquelyAlignedPairs, alignedPairs) << " of aligned)\n";
  cout << "    Multiply : " << multiplyAlignedPairs
       <<"  (" << Fraction(multiplyAlignedPairs, alignedPairs) << " of aligned)\n";
  cout << "  Unaligned  : " << n_pairs - alignedPairs
       <<"  (" << Fraction(n_pairs - alignedPairs, n_pairs) << " of corrected)\n";
  cout << "  Coverage   : " << pairCoverage/genomeSize << "x\n";

  CoverageStatistics(covered, pairCoverage/genomeSize);

  cout << endl;
}
