///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/// ============================================================================

/*

Warning! There is at present a bug in QueryLookupTable, which causes some
alignments to be missed.

This bug arises because for each query sequence, the program finds all
hits to the target sequences (concatenated by MakeLookupTable), and
then clusters these hits by implied start position of the query on the
target.  Unfortunately, if the implied start position is off the
beginning of the target it is interpreted as belonging to the previous
target.  The clustering code should be rewritten to handle such cases
correctly; see the new FindHits code in LookupTable.h, which is
intended to do that, and UniqueUngappedLookup / UngappedLookup /
ShortQueryLookup, which wrap it.

You can reproduce the bug by using the following fasta file, in which
the 3 sequences have identical tails starting at TACCTCTCC...

======= Fasta to reproduce the bug =================================
>Tf7MMP7A
ATGCTACCTCTCCGCGTAGGCGCTCGTTGGTCCAGCAGAGGCGGCCGCCCTTGCGCGAGCAGAA
TGGCGGTAGGGGGTCTAGCTGCGTCTCGTCCGGGGGCATGACACGCAACAGGGGATAGGGACAC
GCACGCAACAGATGG
>Tf9MMP7A
ATGCCGCAAAAACGCAAAACGCAAACGCAACGCATACCTCTCCGCGTAGGCGCTCGTTGGTCCA
GCAGAGGCGGCCGCCCTTGCGCGAGCAGAATGGCGGTAGGGGGTCTAGCTGCGTCTCGTCCGGG
GGCATGACACGCAACAGGGGATAGGGACACGCACGCAACAGATGG
>Tf10MMP7A-Corrected
ATGCATCCTATCCCTATCCCCTATCCCCCTATTACCTCTCCGCGTAGGCGCTCGTTGGTCCAGC
AGAGGCGGCCGCCCTTGCACGAGCAGAATGGCGGTAGGGGGTCTAGCTGCGTCTCGTCCGGGGG
CATGACACGCAACAGGGGATAGGGACACGCACGCAACAGATGG
====================================================================

With default parameters we observe that the 2 vs 2 alignment is not
found, nor are the 1 vs 0 and 1 vs 2 alignments:

--------------------------------------------------------------------------------
Fri Jan 12 14:25:42 2007 run, based on Tue Jan  9 02:22:27 EST 2007 make
QueryLookupTable K=12 SEQS=qlt_bug.fastb L=qlt_bug.lookup
--------------------------------------------------------------------------------
[...]
0fw vs 0, 0 mismatches/0 indels (of 143), from 0-143 to 0-143
0fw vs 1, 1 mismatches/1 indels (of 143), from 0-143 to 29-173
0fw vs 2, 4 mismatches/0 indels (of 143), from 0-143 to 28-171

1fw vs 1, 0 mismatches/0 indels (of 173), from 0-173 to 0-173

2fw vs 1, 18 mismatches/2 indels (of 171), from 0-171 to 0-173
2fw vs 0, 4 mismatches/0 indels (of 171), from 28-171 to 0-143

Adding the SMITH_WAT parameter changes what is found.  This demonstrates
that the lookup table contains enough hits to find the alignments, so
they are being lost in subsequent processing:

--------------------------------------------------------------------------------
Fri Jan 12 14:28:43 2007 run, based on Tue Jan  9 02:22:27 EST 2007 make
QueryLookupTable K=12 SEQS=qlt_bug.fastb L=qlt_bug.lookup SMITH_WAT=True
--------------------------------------------------------------------------------
[...]
0fw vs 0, 0 mismatches/0 indels (of 143), from 0-143 to 0-143
0fw vs 1, 1 mismatches/1 indels (of 143), from 0-143 to 29-173
0fw vs 2, 4 mismatches/0 indels (of 143), from 0-143 to 28-171

1fw vs 1, 0 mismatches/0 indels (of 173), from 0-173 to 0-173

2fw vs 2, 0 mismatches/0 indels (of 171), from 0-171 to 0-171

Now all the k vs k alignments are found, but the 2 vs 1 and 2 vs 0 alignments have disappeared.

The alternate approach using ShortQueryLookup yields the following alignments:
--------------------------------------------------------------------------------
Fri Jan 12 14:47:19 2007 run, based on Thu Jan 11 12:11:05 EST 2007 make
ShortQueryLookup SEQS=qlt_bug.fastb REF=qlt_bug.fastb L=qlt_bug.lookup         \
                 O=qlt_bug MAX_ERRS=20 ERR_DIFF=20
--------------------------------------------------------------------------------
0fw vs 0, 0 mismatches/0 indels (of 143), from 0-143 to 0-143
0fw vs 1, 1 mismatches/1 indels (of 143), from 0-143 to 173-317
0fw vs 2, 4 mismatches/0 indels (of 143), from 0-143 to 345-488

1fw vs 1, 0 mismatches/0 indels (of 173), from 0-173 to 144-317

2fw vs 2, 0 mismatches/0 indels (of 171), from 0-171 to 317-488
2fw vs 1, 18 mismatches/2 indels (of 171), from 0-171 to 144-317

The only alignments now missed are the non-global ones, which is the designed behavior.

*/

/// ============================================================================

/// Find sequences using a k-mer lookup table made by MakeLookupTable.
/// \file QueryLookupTable.cc
///
/// QueryLookupTable.  Find given query sequences in the genome using a
/// k-mer lookup table generated by MakeLookupTable.
///
/// In outline, this works in the following way.
///
/// 1. For each query sequence, find all its locations in the genome which occur
/// with multiplicity at most MAX_FREQ.
///
/// 2. Order these hits by offset, and then break into clusters at places where
/// the offset changes by more than MAX_OFFSET_DIFF.
///
/// 3. Filter these clusters based on size, in a somewhat complex way, determined
/// by parameters MIN_COVERAGE and WINNING_EDGE, yielding "passing" clusters.
///
/// 4. For each offset occurring in a passing cluster, we take all the
/// corresponding hits and try to fill in between them.
///
/// 5. Turn the entire cluster into alignments.
///
/// INPUT PARAMETERS:
///
///    LOOKUP_TABLE (or L): name of lookup table file
///
///    SEQS: name of input file containing query sequence bases
///    (if the name ends with fastb, this is treated as an Arachne-style fastb
///    file; otherwise it must be in fasta format)
///
///    SEQS_IS_FASTB: if True, assume SEQS is a fastb file
///
///    SEQ_NAMES: optional file of names for the query sequences
///
///    QUALS: name of an optional input file containing query sequence quality
///    scores (if the name ends with qualb, this is treated as an Arachne-style
///    qualb file; otherwise it must be in fasta format)
///
///    QUALS_G: quality scores for the genome, as above, optional
///
///    SEQS_TO_PROCESS: if provided, this specifies a subset of the query sequences
///    which are to be processed.  The following formats are allowed:
///    -- @fn, where fn is the name of a file containing some numbers
///    -- n, where n is a number
///    -- {a,b,c} [etc.], where a,b,c are numbers
///    -- [a,b], representing all x with a <= x <= b.
///    The last two formats need to be double-quoted to escape the shell.
///    -- random:n, representing n random read ids (only allowed if SEQS is fastb).
///
///    TARGETS_TO_PROCESS: if provided, this specifies a subset of the target
///    sequences which are to be aligned to.  The format is the same as for
///    SEQS_TO_PROCESS.
///
///    SEQS_TARGETS_TO_PROCESS: if provided, this is the name of a file, whose
///    entries are pairs consisting of a query id and a target id.  (All entries
///    in this file are white space delimited.)  This option cannot be used in
///    conjunction with either SEQS_TO_PROCESS or TARGETS_TO_PROCESS.
///
///    SEQS_TO_EXCLUDE: if provided, this specifies some numeric query sequence ids
///    which are to be excluded (same format as SEQS_TO_PROCESS).
///
/// PERFORMANCE PARAMETERS:
///
///    TMP_DIR: specify directory to be used for writing temporary files
///    (Default: . .)
///
/// HEURISTIC PARAMETERS:
///
///    K: the k-mer size.  This value must be <= the value used to create the
///    lookup table.  Note that using a smaller value than the lookup table value
///    will introduce inefficiencies.
///
///    MAX_FREQ (or MF): k-mers occurring more than this often in the genome are
///    ignored.  More than one value may be specified, separated by colons, in
///    which case the given values are applied in turn, on successively smaller
///    sets of input sequences.  (Default: 500.)
///
///    MAX_OFFSET_DIFF (or MO): for a cluster of k-mer hits to qualify, it must
///    locally have at most this much difference between its offsets.
///    (Default: 1000.)
///
///    MIN_COVERAGE (or MC): fraction of k-mer hits in a cluster must be at least
///    this.  (Default: 0.3.)
///
///    MIN_HITS_TO_OVERRIDE: minimum absolute number of hits in a cluster to
///    override MIN_COVERAGE.  (Default: undefined.)
///
///    WINNING_EDGE (or WE): all clusters are shown that are this close to the best
///    hit, as a multiple of the square root of number of k-mer hits for the best
///    hit.  (Default: 2.5.)
///
///    PROGRESSION_RATIO (or PR): for a cluster of k-mer hits to qualify, the
///    k-mers must cover at least this much of the genome overlap prior to any gap
///    extensions.  (Default: 0.4.)
///
///    KEEP_BEST (or KB): number of highest ranking alignments that will be kept
///    and printed.  (Default: 10.)
///
///    MIN_MUTMER_LENGTH (or MM): minimum length of an extended k-mer hit with
///    errors ("mutmer") for it to be used in generating alignments.
///    (Default: 30.)
///
///    MIN_OVERLAP: minimum length of an overlap between query and target,
///    for it to be accepted regardless of the length of the query sequence.
///    Any overlap covering 20% of the query sequence is automatically accepted.
///    By default, MIN_OVERLAP is undefined.
///
///    SINGLETON_HITS (or SH): if set to True, use a slightly different heuristic
///    for keeping track of k-mer hits.  Observed to be 70% faster sometimes,
///    much slower other times.  (Default: False.)
///
///    END_STRETCH: will attempt to extend alignments to end over k times this
///    number of bases.  (Default: 10.)
///
///    MIN_MATCHES_PERCENT: if specified, reject alignments having less than
///    this percent matching bases over the length of the alignment. Ignores
///    indels, looks only at mismatches.
///
///    MAX_ERROR_PERCENT: if specified, reject alignments having greater than
///    this percent error over the length of the alignment. Ignores
///    indels, looks only at mismatches.
///
///    MIN_BASES_COVERED: if specified, reject alignments having less than
///    this many correct bases in common.
///
///    HEURISTICS: if provided, a file, having one line per iteration, allowing
///    the program to be run multiple times on the same data.  Each line is
///    a list of heuristic parameter definitions as above.  Blank lines and
///    Lines beginning with "#" are ignored.  This option cannot be used in
///    conjunction with the compound form of the MAX_FREQ option.
///
///    SMITH_WAT: if set to True, use banded Smith-Waterman to align.
///
///    BW_ADD: add this value to initial bandwidth for SMITH_WAT option.
///
/// PARSABLE OUTPUT PARAMETERS:
///
///    OUTPUT: full file path where alignments will be printed.  if empty, output
///            will be sent to stdout.
///
///    One or more of the following may be set to True to cause alignments to be
///    printed in the given style.  The default is to have only READABLE_BRIEF
///    set to True.  See README for descriptions.
///
///    PARSEABLE
///
///    PARSEABLE_BRIEF
///
///    VISUAL
///    (note VISUAL_ABBR=True modifies: perfectly aligning stretches of bases
///    are then shown in full)
///
///    READABLE_BRIEF
///
///    RMR_BY_BLOCK
///
///    Naming of target genome contigs (fasta records) is governed by the parameter
///    TARGET_NAMING, which may be set to one of the following values:
///
///         numeric (default): the numeric index of the fasta record, according
///                            to the order encountered by MakeLookupTable
///
///         from_file: take the part of the file name (from which the fasta record
///                    came), after the last slash, and take the part of that
///                    before the first dot, and add [n] to it, where n is the
///                    index of the fasta record in the file
///
///         from_record: the part of the fasta record header which comes after
///         the ">"
///
///         from_record_short: the part of the fasta record header which starts after
///         the ">", and stops at the first white space
///
///         from_record_quoted: the part of the fasta record header which comes
///         after the ">", stripped of double quotes, then double quoted
///
///    Naming of query contigs (fasta records) is governed by the parameter
///    QUERY_NAMING, which may be set to one of the following values:
///
///         numeric (default): the numeric index of the fasta record
///
///         from_record: the part of the read name specified
///         by the command line arg 'nameParser'.  By default, nameParser="first"
///         which grabs the name after the ">" and before the first white space.
///         This reproduces the previous behavior.  If nameParser="last", the name
///         is taken from the last white space until the end of the line.
///
///         from_names_file: from a file specified by the SEQ_NAMES argument
///
///     nameParser, used only with QUERY_NAMING=from_record. Allows different ways
///         to parse a fasta file to get a read name.  See above in 'from_record'.
///
/// LOGGING PARAMETERS:
///
///    PRINT_NQS: if set to True, print statistics on NQS(30,25) discrepancies
///    (requires quality scores for queries and assumes target is high quality)
///
///    PRINT_RMR: if set to True, print "reciprocal match rate (rmr)" for
///    alignment
///
///    PRINT_MM: if set to True, print "match-mismatch score (mm)" for
///    alignment, defined to be (matches-mismatches)/query_length (rounded up to
///    zero if negative)
///
///    PRINT_QUAL_SCORE: if set to True, print Arachne quality score for alignment
///    (you must specify QUALS and QUALS_G to use this)
///
///    MAX_INDELS: if specified, reject alignments having greater than this number
///    of indels
///
///    MAX_MISMATCHES: if specified, reject alignments having greater than this
///    number of mismatches
///
///    MAX_NQS_PERCENT: if specified, reject alignments having greater than this
///    fraction of NQS(30,25) discrepancies
///
///    MAX_QUAL_SCORE: if specified, reject alignments having greater than this
///    Arachne Quality score.  (you must specify QUALS)
///
///    MAX_RMR_PERCENT: if specified, reject alignments having greater than this
///    reciprocal match rate
///
///    MIN_MM_PERCENT: if specified, reject alignments having less than this
///    match-mismatch score
///
///    REQUIRE_PROPER: if specified, require that alignment is proper
///
///    REQUIRE_FULL1: if specified, require that query sequence is aligned from
///    end to end
///
///    REQUIRE_FULL2: if specified, require that target sequence is aligned from
///    end to end (actually only from the second base to the next to the last base,
///    as there seems to be a bug going all the way)
///
///    REQUIRE_POS2_ZERO: if specified, require that alignment start at first base
///    on reference sequence.
///
///    REMOVE_DOMINATED: if set to False, don't filter to remove alignments which
///    dominate others
///
///    MIN_TO_PRINT: if set, only print alignments of queries having at least this
///    many alignments
///
///    FILTER: if set to False, don't filter alignments based on quality
///
///    GLOBAL_STATS: produce global statistics for mismatches and indels
///
///    DUMP_WINDOW: if positive, for each alignment, print the query sequence
///    and the corresponding window on the target sequence, extending by
///    the given amount in both directions.  (This facilitates experimentation
///    with alternate aligners.)
///
///    LIST_UNPLACED: print numerical identifiers of totally unplaced sequences
///
///    UNPLACED_FILE: same as LIST_UNPLACED, but output to given file
///
///    UNPLACED_SEQUENCE_FILE: put unplaced sequences and parts of length >= K
///    in this file
///
///    LIST_UNPLACED_BY_PASS: print totally unplaced sequences after each pass
///
///    ANNOUNCE_ITERATIONS (or AI): announce as each chunk of target genome is
///    processed
///
///    SHOW_PREMUTMERS: if True, show the data from which mutmers are built
///
///    SHOW_MUTMERS: if True, show the mutmers from which alignments are built
///
///    SW_GAP_VERBOSE: if True, turn on verbose logging in the sw_gap aligner
///
///    PRINT_MEMORY_USAGE: if True, print memory usage and related info at each
///    iteration (if ANNOUNCE_ITERATIONS set) and at each major pass.
///
///    NO_HEADER (or NH): if True, omit the initial header showing the name of
///    this executable, etc.
///
///    QUIET: if True, don't print various messages.
///
/// DEBUGGING OPTIONS:
///
///    CHUNKS: if set to a given list, process only the specified chunks of the
///    target genome.  The formatting is the same as for SEQS_TO_PROCESS.
///
///    TRACEBACK_ON_INTERRUPT: if True, provide a traceback if the process is
///    interrupted by a ctrl-C
///
/// OTHER:
///
///    IMPERFECT_EXTENSION: if True, try to extend perfect matches in both
///    directions.
///
///    MAX_PLACEMENTS: ignore query sequences which have more than this number of
///    placements.  The counting of placements is not very intelligent, e.g.
///    if a single query sequence aligns in two nearby parts, it counts as two
///    placements.  Also although this filters the outputting of alignments, it
///    does not prevent the program from "thinking" it has aligned a given query
///    sequence.
///
///    ALIGN_UNALIGNED_BITS: if a query sequence aligns only partially during an
///    iteration, attempt to align the unaligned bits on subsequent iterations.
///    This can be very slow, in part because the bits may be small (and hence
///    harder to align), and in part because even if all iterations have the same
///    parameters, the code is stupid and will keep trying to realign the same bit
///    over and over.
///
///    MIN_QUERY_LENGTH: ignore queries shorter than this.
///
///    SYNC_ALIGNS_TO_TACG: adjust alignments to synchronize with the
///    454-cycles of TACG.
///
///    FILTER_ADD: used in final filter
///
///    FW_ONLY: only allow forward alignments
///    RC_ONLY: only allow reverse alignments
///
///    TRUNCATE_TO: truncate queries to this length
///    TRUNCATE_TO_TAIL: truncate queries to this length, keeping bases at end
///
///    PRINT_PREDICTED_ERRORS: print number of errors predicted by qual scores
///
/// PERFORMANCE NOTES:
///
///    This code can spend a lot of its time reading in the lookup table, over and
///    over.  Moreover, the amount of time spent in doing this is highly dependent
///    on the "state" of the system, including whether there are other i/o
///    intensive processes running and whether this process has been run recently
///    (resulting in possible disk cache speedup), and how many times this
///    process has been run recently.
///
///
/// =================================================================================


#include "Basevector.h"
#include "FastaFilestream.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "ParseSet.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "lookup/FlowAlignSummary.h"
#include "lookup/LookAlign.h"
#include "lookup/LookupTable.h"
#include "lookup/QueryLookupTableCore.h"
#include "math/Arith.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "math/MapIntInt.h"
#include "pairwise_aligners/MakeAlignsMethod.h"
#include "pairwise_aligners/Mutmer.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "random/Random.h"
#include "system/ParsedArgs.h"
#include "system/file/FileReader.h"

#define ABORT(MSG)                                 \
{    out << MSG << "  Abort." << endl << endl;     \
     exit(1);    }

/// Class umutmer: same at mutmer, but the coordinates are unsigned ints.

class umutmer {

     public:

     umutmer( ) { }

     umutmer(unsigned int pos1, unsigned int pos2, int len, int e) :
          pos1_(pos1), pos2_(pos2), len_(len), e_(e) { }

     unsigned int Pos1( ) const { return pos1_; }
     unsigned int Pos2( ) const { return pos2_; }
     int Length( ) const { return len_; }
     int Errors( ) const { return e_; }

     private:

     unsigned int pos1_, pos2_;
     int len_, e_;

};

/// The maximum query size is 200,000,000.  The maximum number of bases
/// in the target fasta files is 4,000,000,000.  These two numbers add up
/// to less than 2^32, for a reason!

const unsigned int max_query_size = 200 * 1000 * 1000;

class index_seq_pos {

     public:

     unsigned int index;
     unsigned int seq;
     unsigned int pos;

     index_seq_pos( ) { }
     index_seq_pos( unsigned int index_arg,
		    unsigned int seq_arg,
		    unsigned int pos_arg )
       : index(index_arg),
	 seq(seq_arg),
	 pos(pos_arg)
     { }

     friend Bool operator<( const index_seq_pos& x1, const index_seq_pos& x2 )
     {    if ( x1.index < x2.index ) return True;
          if ( x1.index > x2.index ) return False;
          if ( x1.seq < x2.seq ) return True;
          return False;    }

};

void FetchHits( unsigned int min_offset, unsigned int max_offset,
     const basevector& s, unsigned int K, unsigned int mfreq,
     vec< pair<unsigned int, unsigned int> >& hitsp, lookup_table& look,
     unsigned int four_to_Kdiff, Bool strict )
{
     hitsp.clear( );
     unsigned int min_on_genome = min_offset;
     unsigned int max_on_genome = max_offset + s.size( ) - K;
     if ( min_on_genome < max_query_size ) min_on_genome = 0;
     else min_on_genome -= max_query_size;
     unsigned int bases_start = look.BasesStart( );
     min_on_genome = Max( bases_start, min_on_genome );
     max_on_genome -= max_query_size;
     max_on_genome
          = Min( max_on_genome, look.BasesStart( ) + look.Bases( ).size( ) - K );
     static vec<index_seq_pos> isp;
     isp.clear( );
     for ( unsigned int u = 0; u <= s.size( ) - K; u++ )
     {    unsigned int index = Index( s, u, K );
          unsigned int freq = 0, index2 = index;
          for ( unsigned int v = 0; v < four_to_Kdiff; v++ )
               freq += look.Freq(index2++);
          if ( freq > mfreq ) continue;
          isp.push_back( index_seq_pos( index, 0, u ) );    }
     for ( unsigned int u = min_on_genome; u <= max_on_genome; u++ )
     {    unsigned int index = Index( look.Bases( ), u - bases_start, K );
          unsigned int freq = 0, index2 = index;
          for ( unsigned int v = 0; v < four_to_Kdiff; v++ )
               freq += look.Freq(index2++);
          if ( freq > mfreq ) continue;
          isp.push_back( index_seq_pos(index, 1, u) );    }
     Sort(isp);
     for ( int u = 0; u < isp.isize( ); u++ )
     {    unsigned int v, w;
          for ( v = u + 1; v < isp.size( ); v++ )
               if ( isp[v].index != isp[u].index ) break;
          if ( isp[u].seq == 1 )
          {    u = v - 1;
               continue;    }
          for ( w = u + 1; w < v; w++ )
               if ( isp[w].seq != isp[u].seq ) break;
          for ( unsigned int x = u; x < w; x++ )
          {    for ( unsigned int y = w; y < v; y++ )
               {    unsigned int offset = isp[y].pos
                         + ( max_query_size - isp[x].pos );
                    if (strict)
                    {    if ( offset < min_offset ) continue;
                         if ( offset > max_offset ) continue;    }
                    hitsp.push_back( make_pair( offset, isp[x].pos ) );    }    }
          u = v - 1;    }
     Sort(hitsp);    }

void ProcessCluster( int call, unsigned int K, int seq_id, Bool rc_seq, int h1,
     int h2, const vec< pair<unsigned int, unsigned int> >& hits, lookup_table& look,
     const basevector& s, int start_on_query, int stop_on_query,
     int orig_query_length, vec<look_align>& qualifiers, int actual_hits,
     Bool SHOW_PREMUTMERS, Bool SHOW_MUTMERS, Bool& have_spoken,
     const map_longlong_longlong& to_query_id, double PROGRESSION_RATIO,
     unsigned int MIN_MUTMER_LENGTH, unsigned int MAX_OFFSET_DIFF,
     unsigned int MIN_OVERLAP, unsigned int END_STRETCH, unsigned int mfreq,
     unsigned int four_to_Kdiff, Bool imperfect_extension,
     const Bool sync_aligns_to_TACG, Bool SMITH_WAT, int BW_ADD,
     const int SW_MISMATCH_PENALTY, const int SW_GAP_PENALTY,
     Bool SW_GAP_VERBOSE, ostream& out )
{
     // Say what we're doing.

     if ( (SHOW_PREMUTMERS || SHOW_MUTMERS) && !have_spoken )
     {    out << "\n*** looking at alignment data for sequence "
               << to_query_id(seq_id) << " ***\n";
          have_spoken = True;    }

     // Determine if all the hits land on the same target contig.  In rare
     // instances, this is not the case, and we have to split up the hits, then
     // call ProcessCluster on each subset.

     const unsigned int UNDEFINED_CONTIG = 2000000000u;
     static unsigned int current_contig(UNDEFINED_CONTIG);
     {    unsigned int start
               = ( hits[h1].first + hits[h1].second ) - max_query_size;
          if ( current_contig == UNDEFINED_CONTIG
               || start < look.ContigStart(current_contig)
               || start >= look.ContigStop(current_contig) )
          {    unsigned int cpos;
               look.GetContigPos( start, current_contig, cpos );    }
          Bool off_contig = False;
          for ( int h = h1; h < h2; h++ )
          {    unsigned int startx
                    = ( hits[h].first + hits[h].second ) - max_query_size;
               if ( startx < look.ContigStart(current_contig)
                    || startx + K > look.ContigStop(current_contig) )
               {    off_contig = True;    }    }
          if (off_contig)
          {    static vec< vec< pair<unsigned int, unsigned int> > > vhits;
               vhits.clear( );
               static vec<int> contigs;
               contigs.clear( );
               static vec<Bool> used;
               used.resize_and_set( h2 - h1, False );
               while(1)
               {    int h0;
                    for ( h0 = h1; h0 < h2; h0++ )
                         if ( !used[ h0 - h1 ] ) break;
                    if ( h0 == h2 ) break;
                    unsigned int cpos;
                    unsigned int startx
                         = ( hits[h0].first + hits[h0].second ) - max_query_size;
                    look.GetContigPos( startx, current_contig, cpos );
                    static vec< pair<unsigned int, unsigned int> > thits;
                    thits.clear( );
                    for ( int h = h0; h < h2; h++ )
                    {    if ( used[ h - h1 ] ) continue;
                         startx = hits[h].first + hits[h].second - max_query_size;
                         if ( startx >= look.ContigStart(current_contig)
                              && startx < look.ContigStop(current_contig) )
                         {    if ( startx + K <= look.ContigStop(current_contig) )
                                   thits.push_back( hits[h] );
                              used[ h - h1 ] = True;    }    }
                    if ( thits.size( ) == 0 ) continue;
                    contigs.push_back(current_contig);
                    vhits.push_back(thits);    }
               for ( int i = 0; i < vhits.isize( ); i++ )
               {    ProcessCluster( call, K, seq_id, rc_seq, 0, vhits[i].size( ),
                         vhits[i], look, s, start_on_query, stop_on_query,
                         orig_query_length, qualifiers, actual_hits,
                         SHOW_PREMUTMERS, SHOW_MUTMERS, have_spoken, to_query_id,
                         PROGRESSION_RATIO, MIN_MUTMER_LENGTH, MAX_OFFSET_DIFF,
                         MIN_OVERLAP, END_STRETCH, mfreq, four_to_Kdiff,
                         imperfect_extension, sync_aligns_to_TACG, SMITH_WAT,
                         BW_ADD, SW_MISMATCH_PENALTY, SW_GAP_PENALTY, 
                         SW_GAP_VERBOSE, out );    }
               return;    }    }

     // Create merged hits, which are presented as umutmers.

     static vec<umutmer> um;
     um.clear( );
     for ( int h = h1; h < h2; h++ )
     {    int x;
          for ( x = h + 1; x < h2; x++ )
               if ( hits[x].first != hits[h].first ) break;

          // The hits between h and x all have the same offset.

          unsigned int offset = hits[h].first;
          unsigned int start = hits[h].second, stop = hits[x-1].second + K;

          // The bases on the sequence from start to stop all have the potential to
          // match at the given offset.  First determine which ones actually match.

          static vec<Bool> match;
          match.resize_and_set( s.size( ), False );
          for ( unsigned int y = start; y < stop; y++ )
          {    if ( s[y] == look.Base( (offset + y) - max_query_size ) )
                    match[y] = True;    }

          // Do perfect extensions off the ends.

          Bool off_contig_end = False;
          if ( start != 0 )
          {    for ( unsigned int y = start - 1; ; y-- )
               {    if ( offset + y < max_query_size ) break;
                    unsigned int bl = (offset + y) - max_query_size;
                    if ( !look.BaseInMemory(bl) ||
                         bl < look.ContigStart(current_contig) )
                    {    off_contig_end = True;
                         break;    }
                    if ( s[y] != look.Base(bl) ) break;
                    match[y] = True;
                    --start;
                    if ( y == 0 ) break;    }    }
          for ( unsigned int y = stop; y < s.size( ); y++ )
          {    unsigned int bl = (offset + y) - max_query_size;
               if ( !look.BaseInMemory(bl) || bl >= look.ContigStop(current_contig) )
               {    off_contig_end = True;
                    break;    }
               if ( s[y] != look.Base(bl) ) break;
               match[y] = True;
               ++stop;    }

          // Do imperfect extensions off the ends.

          unsigned int start0 = start, stop0 = stop;
          int imperfect_left_extension = 0, imperfect_right_extension = 0;
          int imperfect_left_extension_errs = 0, imperfect_right_extension_errs = 0;
          if (imperfect_extension)
          {    const int peek_floor = -2;
               Bool progress;
               do
               {    progress = False;
                    static vec<int> peek;
                    peek.clear( );
                    int peeksum = 0;
                    if ( start != 0 )
                    {    for ( unsigned int y = start - 1; ; y-- )
                         {    if ( offset + y < max_query_size ) break;
                              unsigned int bl = (offset + y) - max_query_size;
                              if ( !look.BaseInMemory(bl) ||
                                   bl < look.ContigStart(current_contig) )
                              {    break;    }
                              if ( s[y] == look.Base(bl) ) peek.push_back(+1);
                              else peek.push_back(-1);
                              peeksum += peek.back( );
                              if ( peeksum < peek_floor ) break;
                              if ( peeksum > 0 ) break;
                              if ( y == 0 ) break;    }    }
                    if ( peeksum > 0 )
                    {    progress = True;
                         unsigned int count = 0;
                         for ( unsigned int y = start - 1; ; y-- )
                         {    match[y] = ( peek[count] > 0 );
                              if ( peek[count] < 0 )
                                   imperfect_left_extension_errs++;
                              if ( ++count == peek.size( ) ) break;    }
                         imperfect_left_extension += peek.size( );
                         start -= (int) peek.size( );    }
                    peek.clear( );
                    peeksum = 0;
                    for ( unsigned int y = stop; y < s.size( ); y++ )
                    {    unsigned int bl = (offset + y) - max_query_size;
                         if ( !look.BaseInMemory(bl)
                              || bl >= look.ContigStop(current_contig) )
                         {    break;    }
                         if ( s[y] == look.Base(bl) ) peek.push_back(+1);
                         else peek.push_back(-1);
                         peeksum += peek.back( );
                         if ( peeksum < peek_floor ) break;
                         if ( peeksum > 0 ) break;    }
                    if ( peeksum > 0 )
                    {    progress = True;
                         unsigned int count = 0;
                         for ( unsigned int y = stop; y < s.size( ); y++ )
                         {    match[y] = ( peek[count] > 0 );
                              if ( peek[count] < 0 )
                                   imperfect_right_extension_errs++;
                              if ( ++count == peek.size( ) ) break;    }
                         imperfect_right_extension += peek.size( );
                         stop += (int) peek.size( );    }    }
               while(progress);    }

          // Find the perfect blocks of length >= K.  Stored as {(start, length)}.

          static vec< pair<unsigned int, unsigned int> > perfect;
          perfect.clear( );
          unsigned int last = start0;
          for ( unsigned int y = start0; y < stop0; y++ )
          {    if ( !match[y] )
               {    if ( y - last >= K )
                         perfect.push_back( make_pair( last, y - last ) );
                    last = y + 1;    }    }
          if ( stop0 - last >= K )
               perfect.push_back( make_pair(last, stop0 - last) );

          // Compute length and number of errors for each gap
          // between the length >= K perfect blocks.

          static vec< pair<int, int> > between;
          between.resize( perfect.size( ) - 1 );
          for ( int u = 0; u < between.isize( ); u++ )
          {    int length = perfect[u+1].first
                    - perfect[u].first - perfect[u].second;
               int errors = 0, st = perfect[u].first + perfect[u].second;
               for ( int v = 0; v < length; v++ )
               {    if ( s[ st + v ]
                         != look.Base( (offset + st + v) - max_query_size ) )
                    {    ++errors;    }    }
               between[u] = make_pair( length, errors );    }

          if (SHOW_PREMUTMERS)
          {    out << "\nstart=" << perfect[0].first << ": ";
               for ( int u = 0; u < perfect.isize( ); u++ )
               {    out << "perf=" << perfect[u].second;
                    if ( u < (int) perfect.size( ) - 1 )
                         out << "/" << between[u].second << "_errs_of_"
                              << between[u].first << "/";    }
               out << "\n";    }

          // Decide to close gap if it has <= 5 errors or <= 20% error rate.

          static vec<Bool> to_close;
          to_close.resize( between.size( ) );
          for ( int u = 0; u < between.isize( ); u++ )
               to_close[u] = between[u].second <= 5
                    || float( between[u].second ) / float( between[u].first ) <= 0.2;

          // Close gaps, yielding umutmers.

          last = 0;
          for ( int u = 0; u <= to_close.isize( ); u++ )
          {    if (u == (int) to_close.size( ) || !to_close[u])
               {    unsigned int pos1 = perfect[last].first;
                    int len = int( ( perfect[u].first + perfect[u].second ) - pos1 );
                    unsigned int pos2 = (pos1 + offset) - max_query_size;
                    int e = 0;
                    for ( int v = last; v < u; v++ )
                         e += between[v].second;
                    if ( u == 0 )
                    {    pos1 -= imperfect_left_extension;
                         pos2 -= imperfect_left_extension;
                         len += imperfect_left_extension;
                         e += imperfect_left_extension_errs;    }
                    if ( u == (int) to_close.size( ) )
                    {    len += imperfect_right_extension;
                         e += imperfect_right_extension_errs;    }
                    um.push_back( umutmer( pos1, pos2, len, e ) );
                    if ( len >= (int) MIN_MUTMER_LENGTH )
                    {    if (SHOW_MUTMERS)
                         {    out << "[offset "
                                   << (longlong) pos2 - (longlong) pos1
                                   << ", hits = " << x - h
                                   << ", pos1 = " << pos1 << ", len = " << len
                                   << ", errs = " << e << "]\n";    }    }
                    last = u + 1;    }    }

          h = x - 1;    }

     // Convert umutmers to mutmers, try to build alignment.  Note that offset
     // values we compute now are (as always) shifted by max_query_size to prevent
     // them from going negative inside an unsigned int.

     if ( um.size( ) == 0 ) return;
     unsigned int min_offset = ( max_query_size + um[0].Pos2( ) ) - um[0].Pos1( );
     unsigned int max_offset = ( max_query_size + um[0].Pos2( ) ) - um[0].Pos1( );
     for ( int r = 1; r < um.isize( ); r++ )
     {    min_offset = Min( min_offset,
               ( max_query_size + um[r].Pos2( ) ) - um[r].Pos1( ) );
          max_offset =
               Max( max_offset,
                    ( max_query_size + um[r].Pos2( ) ) - um[r].Pos1( ) );    }

     // (Make sure we're not going off the end of the chunk of bases in memory.)

     unsigned int extra = 100;
     longlong start0 = ( (longlong)(min_offset) - (longlong)(extra) )
          - (longlong) max_query_size;
     longlong stop0 = (longlong)( max_offset + s.size( ) + extra )
          - (longlong) max_query_size;
     unsigned int start
          = Max( start0, (longlong) look.BasesStart( ),
               (longlong) look.ContigStart(current_contig) );
     unsigned int stop
          = Min( stop0, (longlong) look.BasesStop( ),
               (longlong) look.ContigStop(current_contig) );

     // Test for a rare failure mode.  Not sure how this can happen.  Not a very
     // good solution.

     if ( !( stop > start ) ) return;

     ForceAssertGe( start, look.BasesStart( ) );

     static basevector bpart;
     bpart.Setsize( stop - start );
     bpart.SetToSubOf( look.Bases( ), start - look.BasesStart( ), stop - start );

     // Test for a rare failure mode.  Not sure how this can happen.  Not a very
     // good solution.

     Bool neg_um = False;
     for ( int r = 0; r < um.isize( ); r++ )
          if ( um[r].Pos2( ) < start ) neg_um = True;
     if (neg_um) return;

     static vec<mutmer> mu;
     mu.resize( um.size( ) );
     for ( int r = 0; r < um.isize( ); r++ )
          mu[r].SetFrom( um[r].Pos1( ), (int) ( um[r].Pos2( ) - start ),
               um[r].Length( ), um[r].Errors( ) );
     for ( int r = 0; r < mu.isize( ); r++ )
     {    ForceAssertGe( mu[r].Pos1( ), 0 );
          ForceAssertGe( mu[r].Pos2( ), 0 );    }

     // Test for a rare failure mode.  Not sure how this can happen.  Not a very
     // good solution.

     Bool big_um = False;
     for ( int r = 0; r < mu.isize( ); r++ )
          if ( mu[r].Pos2( ) + mu[r].Length( ) > (int) bpart.size( ) ) big_um = True;
     if (big_um) return;

     if (SMITH_WAT)
     {    int offset_low = mu[0].Offset( ), offset_high = mu[0].Offset( );
          for ( int r = 1; r < mu.isize( ); r++ )
          {    offset_low = Min( offset_low, mu[r].Offset( ) );
               offset_high = Max( offset_high, mu[r].Offset( ) );    }
          int offset = (offset_low + offset_high) / 2;
          int bandwidth = (offset_high - offset_low) / 2;
          static align a;
          int errs;
          SmithWatBandedA( s, bpart, offset, bandwidth + BW_ADD, a, errs,
               0, SW_MISMATCH_PENALTY, SW_GAP_PENALTY );
          if ( sync_aligns_to_TACG ) a.Sync_to_TACG( s, bpart, rc_seq );
          static vector<int> mgg;
          mgg = a.MutationsGap1Gap2( s, bpart );
          int mutations = mgg[0], indels = mgg[1] + mgg[2];
          unsigned int start2 = start + (unsigned int) a.pos2( );
          unsigned int c, cpos;
          look.GetContigPos( start2, c, cpos );
          a.Setpos2(cpos);
          if ( !rc_seq ) a.AddToPos1(start_on_query);
          else a.AddToPos1( orig_query_length - stop_on_query );
          qualifiers.push_back( look_align( seq_id, c, orig_query_length,
               look.ContigSize(c), rc_seq, a, actual_hits, mutations, indels ) );
          return;    }


     // Define makealigns method.

     makealigns_method *method_ptr;
     makealigns_sw_gap_method sw_gap_method;
     sw_gap_method.SetMaxErrs(10000);
     sw_gap_method.SetEndStretch(END_STRETCH);
     sw_gap_method.SetMinMaxMutmerLength(MIN_MUTMER_LENGTH);
     sw_gap_method.SetMaxGap(400);
     sw_gap_method.SetMinProgressionLength(0);
     sw_gap_method.SetMinProgressionRatio(PROGRESSION_RATIO);
     sw_gap_method.SetMinOverlapFraction(0.2);
     sw_gap_method.SetIgnoreOverlapFractionLength(MIN_OVERLAP);
     sw_gap_method.SetVerbose(SW_GAP_VERBOSE);
     sw_gap_method.SetAffinePenalties(True);
     method_ptr = &sw_gap_method;

     // Build alignments.

     static vec<align> aligns(1000);
     static vec<int> errors(1000);
     int aligns_length = 0, min_mutmer = Min( 12, (int) MIN_MUTMER_LENGTH );
     static vec<int> od;
     od.clear( );
     od.push_back(MAX_OFFSET_DIFF);
     if ( MAX_OFFSET_DIFF > 100 ) od.push_back(100);
     unsigned int max_pos1 = s.size( ), min_Pos1 = 0;
     Bool found_full = False;
     for ( int mpass = 0; mpass < od.isize( ); mpass++ )
     {    sw_gap_method.SetMaxMutmerOffsetDiff( od[mpass] );
          Bool answer = method_ptr->MutmersToAlign( mu, K, s, bpart, aligns, errors,
               aligns_length, min_mutmer, ( SW_GAP_VERBOSE ? &out : 0 ) );
          if ( !answer ) aligns_length = 0;

          // Find unused mutmers.

          static vec<Bool> mu_used;
          mu_used.resize_and_set( mu.size( ), False );
          for ( int u = 0; u < aligns_length; u++ )
          {    const align& a = aligns[u];
               int p1 = a.pos1( ), p2 = a.pos2( );
               for ( int j = 0; j < a.Nblocks( ); j++ )
               {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
                    if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
                    for ( int k = 0; k < mu.isize( ); k++ )
                    {    if ( mu[k].Pos1( ) - mu[k].Pos2( ) != p1 - p2 ) continue;
                         int over = IntervalOverlap( p1, p1 + a.Lengths(j),
                              mu[k].Pos1( ), mu[k].Pos1( ) + mu[k].Length( ) );
                         if ( float(over)/float(mu[k].Length( )) < 0.75 ) continue;
                         mu_used[k] = True;    }
                    p1 += a.Lengths(j);
                    p2 += a.Lengths(j);    }    }
          static vec<mutmer> mu2;
          mu2.clear( );
          for ( int i = 0; i < mu_used.isize( ); i++ )
               if ( !mu_used[i] ) mu2.push_back( mu[i] );
          mu = mu2;

          // Save alignments.

          for ( int u = 0; u < aligns_length; u++ )
          {    align& a = aligns[u];

	       // adjust the align object if necessary
	       if ( sync_aligns_to_TACG )
		 {
		   a.Sync_to_TACG( s, bpart, rc_seq );
		 }

               static vector<int> mgg;
               mgg = a.MutationsGap1Gap2( s, bpart );
               int mutations = mgg[0], indels = mgg[1] + mgg[2];
               unsigned int start2 = start + (unsigned int) a.pos2( );
               unsigned int c, cpos;
               look.GetContigPos( start2, c, cpos );
               a.Setpos2(cpos);
               if ( !rc_seq ) a.AddToPos1(start_on_query);
               else a.AddToPos1( orig_query_length - stop_on_query );
               qualifiers.push_back( look_align( seq_id, c, orig_query_length,
                    look.ContigSize(c), rc_seq, a, actual_hits, mutations,
                    indels ) );
               if ( a.pos1( ) == 0 && a.Pos1( ) == (int) s.size( ) )
                    found_full = True;
               max_pos1 = Max( max_pos1, (unsigned int) a.pos1( ) );
               min_Pos1 = Min( min_Pos1, (unsigned int) a.Pos1( ) );    }    }

     // If only partial alignments were found, try doing an extended search for
     // hits in the neighborhood, then retrying.
     //
     // Note (7/14/04).  This was EXTREMELY slow in certain situations, so I
     // added the restriction that s.size( ) < 2000.  I'm not sure what the right
     // solution is.

     if ( call == 0 && !found_full && ( max_pos1 > 0 || min_Pos1 < s.size( ) )
          && s.size( ) < 2000 )
     {    int mfreq_mult = 600;
          unsigned int min_offset = hits[h1].first, max_offset = hits[h1].first;
          for ( int h = h1+1; h < h2; h++ )
          {    min_offset = Min( min_offset, hits[h].first );
               max_offset = Max( max_offset, hits[h].first );    }
          if ( min_offset >= max_pos1 ) min_offset -= max_pos1;
          else min_offset = 0;
          max_offset += (int) s.size( ) - min_Pos1;
          static vec< pair<unsigned int, unsigned int> > hits2;
          FetchHits( min_offset, max_offset, s, K, mfreq_mult * mfreq, hits2, look,
               four_to_Kdiff, False );
          ProcessCluster( 1, K, seq_id, rc_seq, 0, hits2.size( ), hits2, look, s,
               start_on_query, stop_on_query, orig_query_length, qualifiers,
               actual_hits, SHOW_PREMUTMERS, SHOW_MUTMERS, have_spoken, to_query_id,
               PROGRESSION_RATIO, MIN_MUTMER_LENGTH, MAX_OFFSET_DIFF,
               MIN_OVERLAP, END_STRETCH, mfreq, four_to_Kdiff, imperfect_extension,
	       sync_aligns_to_TACG, SMITH_WAT, BW_ADD, SW_MISMATCH_PENALTY,
               SW_GAP_PENALTY, SW_GAP_VERBOSE, out );    }    }

void ProcessQuerySequence( const int npasses, const Bool firstpass, 
     lookup_table& look, int seq_id,
     basevector s, int start_on_query, int stop_on_query, int orig_query_length,
     qualvector q, int& n_all_qualifiers, off_t& disk_pos, unsigned int Kdiffbits,
     unsigned int four_to_Kdiff, unsigned int mfreq, const map_longlong_longlong& to_query_id,
     unsigned int K, double MIN_COVERAGE, unsigned int MIN_HITS_TO_OVERRIDE,
     double WINNING_EDGE, unsigned int MAX_OFFSET_DIFF, unsigned int MIN_OVERLAP,
     unsigned int END_STRETCH, double PROGRESSION_RATIO,
     unsigned int MIN_MUTMER_LENGTH, Bool SINGLETON_HITS, unsigned int KEEP_BEST,
     Bool SHOW_PREMUTMERS, Bool SHOW_MUTMERS, vec<unsigned int>* hits_ptr,
     vec< pair<unsigned int, unsigned int> >* hitsp_ptr, int fdq, int fdqs,
     Bool imperfect_extension, const vec<longlong>& targets_to_process,
     const Bool sync_aligns_to_TACG, Bool SMITH_WAT, int BW_ADD,
     const int SW_MISMATCH_PENALTY, const int SW_GAP_PENALTY,
     Bool SW_GAP_VERBOSE, ostream& out )
{
     vec<unsigned int>& hits = *hits_ptr;
     vec< pair<unsigned int, unsigned int> >& hitsp = *hitsp_ptr;
     static vec<look_align> qualifiers;
     qualifiers.clear( );
     int min_hits = int(floor(
          MIN_COVERAGE * ((double) (stop_on_query-start_on_query) - (double) K) ) );
     min_hits = Min( min_hits, (int) MIN_HITS_TO_OVERRIDE );
     int winning_edge = int( floor( WINNING_EDGE * sqrtf(float(min_hits)) ) );
     winning_edge = Min( winning_edge, min_hits );
     int hfloor = min_hits - winning_edge;

     if ( start_on_query > 0 || stop_on_query < (int) s.size( ) )
     {    static basevector t;
          t.SetToSubOf( s, start_on_query, stop_on_query - start_on_query );
          s = t;    }

     for ( int pass = firstpass; pass <= npasses; pass++ )
     {    if ( pass == 2 )
          {    s.ReverseComplement( );
               q.ReverseMe( );    }

          Bool have_spoken = False;

          if ( hits_ptr != 0 ) hits.clear( );
          if ( hitsp_ptr != 0 ) hitsp.clear( );
          int most_hits_in_cluster = 0;

          // Generate hits.  These are pairs (offset, qpos), consisting of offset =
          // the difference between target sequence position and query sequence
          // position, and qpos = query sequence position.  The offset values are
          // shifted by max_query_size, so that we do not attempt to represent
          // negative numbers in an unsigned int.
          //
          // Note.  If SINGLETON_HITS=True, we generate only the offsets now (hits).
          // Then we generate hitsp, whose entries are pairs.

          for ( unsigned int r = 0; r <= s.size( ) - K; r++ )
          {    unsigned int index = Index( s, r, K ) << Kdiffbits;
               unsigned int freq = 0, index2 = index;
               for ( unsigned int v = 0; v < four_to_Kdiff; v++ )
                    freq += look.Freq(index2++);
               if ( freq > mfreq ) continue;
               for ( unsigned int v = 0; v < four_to_Kdiff; v++ )
               {    unsigned int start = look.StartLocs(index);
                    unsigned int stop = look.StopLocs(index);
                    if ( targets_to_process.empty( ) )
                    {    if (SINGLETON_HITS)
                         {    for ( unsigned int l = start; l < stop; l++ )
                              {    unsigned int offset
                                        = look.Locs(l) + ( max_query_size - r );
                                   hits.push_back(offset);    }    }
                         else
                         {    for ( unsigned int l = start; l < stop; l++ )
                              {    unsigned int offset
                                        = look.Locs(l) + ( max_query_size - r );
                                   hitsp.push_back(
                                        make_pair( offset, r ) );    }    }    }
                    else
                    {    if (SINGLETON_HITS)
                         {    for ( unsigned int l = start; l < stop; l++ )
                              {    unsigned int loc = look.Locs(l);
                                   Bool in_targets = False;
                                   for ( int u = 0; u < targets_to_process.isize( );
                                        u++ )
                                   {    int t = targets_to_process[u];
                                        if ( loc >= look.ContigStart(t)
                                             && loc <= look.ContigStop(t) )
                                        {    in_targets = True;    }    }
                                   if ( !in_targets ) continue;
                                   unsigned int offset = loc + (max_query_size - r);
                                   hits.push_back(offset);    }    }
                         else
                         {    for ( unsigned int l = start; l < stop; l++ )
                              {    unsigned int loc = look.Locs(l);
                                   Bool in_targets = False;
                                   for ( int u = 0; u < targets_to_process.isize( );
                                        u++ )
                                   {    int t = targets_to_process[u];
                                        if ( loc >= look.ContigStart(t)
                                             && loc <= look.ContigStop(t) )
                                        {    in_targets = True;    }    }
                                   if ( !in_targets ) continue;
                                   unsigned int offset = loc + (max_query_size - r);
                                   hitsp.push_back(
                                        make_pair( offset, r ) );    }    }    }
                    ++index;    }    }

          if ( hits_ptr != 0 ) Sort(hits);
          if ( hitsp_ptr != 0 ) Sort(hitsp);
          int total_hits = ( SINGLETON_HITS ? hits.size( ) : hitsp.size( ) );
          for ( int h1 = 0; h1 < total_hits; h1++ )
          {    int h2;

               if (SINGLETON_HITS)
               {    for ( h2 = h1 + 1; h2 < hits.isize( ); h2++ )
                    {    if ( hits[h2] - hits[h1] > MAX_OFFSET_DIFF) break;    }    }
               else
               {    for ( h2 = h1 + 1; h2 < hitsp.isize( ); h2++ )
                    {    if ( hitsp[h2].first - hitsp[h1].first
                              > MAX_OFFSET_DIFF ) break;    }    }

               // If the number of hits is below the threshold, ignore.

               if ( h2 - h1 < hfloor )
               {    h1 = h2 - 1;
                    continue;    }

               // Compute hitspl.

               static vec< pair<unsigned int, unsigned int> > hitspl;
               if (SINGLETON_HITS)
               {    int mfreq_mult = 4;
                    unsigned int min_offset = hits[h1], max_offset = hits[h1];
                    for ( int h = h1+1; h < h2; h++ )
                    {    min_offset = Min( min_offset, hits[h] );
                         max_offset = Max( max_offset, hits[h] );    }
                    FetchHits( min_offset, max_offset, s, K, mfreq_mult * mfreq,
                         hitspl, look, four_to_Kdiff, True );    }

               // Compute number of actual hits.


               static vec<Bool> seqhits;
               seqhits.resize_and_set( s.size( ) - K + 1, False );
               if (SINGLETON_HITS)
               {    for ( unsigned int h = 0; h < hitspl.size( ); h++ )
                    {    int u = hitspl[h].second;
                         unsigned int index = Index( s, u, K );
                         unsigned int freq = 0, index2 = index;
                         for ( unsigned int v = 0; v < four_to_Kdiff; v++ )
                              freq += look.Freq(index2++);
                         if ( freq > mfreq ) continue;
                         seqhits[u] = True;    }    }
               else
               {    for ( int h = h1; h < h2; h++ )
                         seqhits[ hitsp[h].second ] = True;    }
               int actual_hits = Sum(seqhits);

               // If the number of actual hits is below the threshold, ignore.

               if ( actual_hits < hfloor )
               {    h1 = h2 - 1;
                    continue;    }

               most_hits_in_cluster = Max( most_hits_in_cluster, actual_hits );
               hfloor = Max( hfloor,
                    most_hits_in_cluster -
                         int( floor(
                         WINNING_EDGE * sqrtf(float(most_hits_in_cluster)) ) ) );
               if ( hfloor < 0 ) hfloor = 0;

               // Generate alignments.

               const vec< pair<unsigned int, unsigned int> >& hp
                    = ( SINGLETON_HITS ? hitspl : hitsp );
               int hp1 = ( SINGLETON_HITS ? 0 : h1 );
               int hp2 = ( SINGLETON_HITS ? (int) hp.size( ) : h2 );
               ProcessCluster( 0, K, seq_id, pass == 2, hp1, hp2, hp, look, s,
                    start_on_query, stop_on_query, orig_query_length, qualifiers,
                    actual_hits, SHOW_PREMUTMERS, SHOW_MUTMERS, have_spoken,
                    to_query_id, PROGRESSION_RATIO, MIN_MUTMER_LENGTH,
                    MAX_OFFSET_DIFF, MIN_OVERLAP, END_STRETCH, mfreq, four_to_Kdiff,
                    imperfect_extension, sync_aligns_to_TACG, SMITH_WAT, BW_ADD,
                    SW_MISMATCH_PENALTY, SW_GAP_PENALTY, SW_GAP_VERBOSE, out );

               h1 = h2 - 1;    }    }

     // Sort and then prune any alignments ranking worse than KEEP_BEST.

     Sort(qualifiers);

     int qualifiers_saved = 0;
     for ( int m = 0; m < qualifiers.isize( ); m++ )
     {    const look_align& q = qualifiers[m];
          if ( q.nhits >= hfloor )
          {
               // Story query id as longlong to make it easier to
               // read the entire structure back into memory.

               longlong qid = q.query_id;
               ForceAssertEq( sizeof(qid), sizeof(disk_pos) );

               WriteBytes( fdqs, &qid, sizeof(qid) );
               WriteBytes( fdqs, &disk_pos, sizeof(disk_pos) );
               disk_pos += q.BinaryWrite(fdq);
               ++n_all_qualifiers;
               if ( ++qualifiers_saved == (int) KEEP_BEST ) break;     }    }    }

void arachne_signal_handler( int signal_number );

const unsigned int undefined = 1000000000;

void QueryLookupTableCore( int argc, char *argv[] )
{
     BeginCommandArguments;
     CommandArgument_String(SEQS);
     CommandArgument_String_OrDefault(SEQ_NAMES, "");
     CommandArgument_String_OrDefault(QUALS, "");
     CommandArgument_String_OrDefault(QUALS_G, "");
     CommandArgument_UnsignedInt(K);
     CommandArgument_String_Abbr_OrDefault(LOOKUP_TABLE, L, "");
     CommandArgument_UnsignedInt_Abbr_OrDefault(MAX_OFFSET_DIFF, MO, 1000);
     CommandArgument_UnsignedInt_OrDefault(MIN_OVERLAP, undefined);
     CommandArgument_String_Abbr_OrDefault(MAX_FREQ, MF, "500");
     CommandArgument_Double_Abbr_OrDefault(WINNING_EDGE, WE, 2.5);
     CommandArgument_Double_Abbr_OrDefault(MIN_COVERAGE, MC, 0.3);
     CommandArgument_UnsignedInt_OrDefault(MIN_HITS_TO_OVERRIDE, undefined);
     CommandArgument_Double_Abbr_OrDefault(PROGRESSION_RATIO, PR, 0.4);
     CommandArgument_Double_OrDefault(MAX_NQS_PERCENT, undefined);
     CommandArgument_Double_OrDefault(MAX_QUAL_SCORE, undefined);
     CommandArgument_Double_OrDefault(MAX_RMR_PERCENT, undefined);
     CommandArgument_Double_OrDefault(MIN_MM_PERCENT, undefined);
     CommandArgument_Int_OrDefault(MIN_PERFECT_MATCH, 0);
     CommandArgument_UnsignedInt_Abbr_OrDefault(KEEP_BEST, KB, 10);
     CommandArgument_UnsignedInt_Abbr_OrDefault(MIN_MUTMER_LENGTH, MM, 30);
     CommandArgument_Bool_Abbr_OrDefault(SINGLETON_HITS, SH, False);
     CommandArgument_UnsignedInt_OrDefault(END_STRETCH, 10);
     CommandArgument_String_OrDefault(HEURISTICS, "");
     CommandArgument_Bool_OrDefault(SHOW_PREMUTMERS, False);
     CommandArgument_Bool_OrDefault(SHOW_MUTMERS, False);
     CommandArgument_Bool_OrDefault(SMITH_WAT, False);
     CommandArgument_UnsignedInt_OrDefault(BW_ADD, 10);
     CommandArgument_Bool_OrDefault(LIST_UNPLACED, False);
     CommandArgument_Bool_OrDefault(LIST_UNPLACED_BY_PASS, False);
     CommandArgument_Bool_OrDefault(REQUIRE_PROPER, False);
     CommandArgument_Bool_OrDefault(REQUIRE_FULL1, False);
     CommandArgument_Bool_OrDefault(REQUIRE_FULL2, False);
     CommandArgument_Bool_OrDefault(REQUIRE_POS2_ZERO, False);
     CommandArgument_Bool_OrDefault(REMOVE_DOMINATED, True);
     CommandArgument_String_OrDefault(UNPLACED_FILE, "");
     CommandArgument_Bool_Abbr_OrDefault(ANNOUNCE_ITERATIONS, AI, False);
     CommandArgument_UnsignedInt_OrDefault(DUMP_WINDOW, 0);
     CommandArgument_String_OrDefault(OUTPUT, "");
     CommandArgument_Bool_OrDefault(PARSEABLE, False);
     CommandArgument_Bool_OrDefault(PARSEABLE_BRIEF, False);
     CommandArgument_Bool_OrDefault(VISUAL, False);
     CommandArgument_Bool_OrDefault(VISUAL_ABBR, True);
     CommandArgument_Bool_OrDefault(READABLE_BRIEF, True);
     CommandArgument_Bool_OrDefault(RMR_BY_BLOCK, False);
     CommandArgument_Bool_OrDefault(GLOBAL_STATS, False);
     CommandArgument_String_OrDefault(TARGET_NAMING, "numeric");
     CommandArgument_String_OrDefault(QUERY_NAMING, "numeric");
     CommandArgument_String_OrDefault(TARGET_FASTA, "");
     CommandArgument_Bool_OrDefault(PRINT_MEMORY_USAGE, False);
     CommandArgument_String_Abbr_OrDefault(SEQS_TO_PROCESS, STP, "undefString");
     CommandArgument_String_OrDefault(TARGETS_TO_PROCESS, "undefString");
     CommandArgument_String_OrDefault(SEQS_TARGETS_TO_PROCESS, "undefString");
     CommandArgument_String_OrDefault(SEQS_TO_EXCLUDE, "undefString");
     CommandArgument_String_OrDefault(TMP_DIR, ".");
     CommandArgument_String_OrDefault(CHUNKS, "undefString");
     CommandArgument_Bool_OrDefault(IMPERFECT_EXTENSION, False);
     CommandArgument_Bool_OrDefault(FILTER, True);
     CommandArgument_Bool_OrDefault(PRINT_NQS, False);
     CommandArgument_Bool_OrDefault(PRINT_PREDICTED_ERRORS, False);
     CommandArgument_Bool_OrDefault(PRINT_RMR, False);
     CommandArgument_Bool_OrDefault(PRINT_MM, False);
     CommandArgument_Bool_OrDefault(PRINT_QUAL_SCORE, False);
     CommandArgument_Bool_OrDefault(SW_GAP_VERBOSE, False);
     CommandArgument_Int_OrDefault(SW_MISMATCH_PENALTY, 2);
     CommandArgument_Int_OrDefault(SW_GAP_PENALTY, 3);
     CommandArgument_Bool_OrDefault(TRACEBACK_ON_INTERRUPT, False);
     CommandArgument_UnsignedInt_OrDefault(MAX_PLACEMENTS, undefined);
     CommandArgument_String_OrDefault(nameParser, "first");
     CommandArgument_Bool_OrDefault(ALIGN_UNALIGNED_BITS, False);
     CommandArgument_String_OrDefault(UNPLACED_SEQUENCE_FILE, "");
     CommandArgument_UnsignedInt_OrDefault(MAX_MISMATCHES, undefined);
     CommandArgument_UnsignedInt_OrDefault(MAX_INDELS, undefined);
     CommandArgument_Int_OrDefault(MIN_MATCHES_PERCENT, 0);
     CommandArgument_Int_OrDefault(MAX_ERROR_PERCENT, 100);
     CommandArgument_Int_OrDefault(MIN_BASES_COVERED, 0);
     CommandArgument_UnsignedInt_OrDefault(MIN_QUERY_LENGTH, 0);
     CommandArgument_Bool_OrDefault(SYNC_ALIGNS_TO_TACG, False);
     CommandArgument_Bool_OrDefault(QUIET, False);
     CommandArgument_Bool_OrDefault(SEQS_IS_FASTB, False);
     CommandArgument_Bool_OrDefault(QUALS_IS_QUALB, False);
     CommandArgument_Int_OrDefault(FILTER_ADD, 8);
     CommandArgument_Int_OrDefault(FILTER_MULT, 2);
     CommandArgument_String_OrDefault(OUTFILE, "");
     CommandArgument_Bool_OrDefault(FW_ONLY, False);
     CommandArgument_Bool_OrDefault(RC_ONLY, False);
     CommandArgument_Int_OrDefault(TRUNCATE_TO, -1);
     CommandArgument_Int_OrDefault(TRUNCATE_TO_TAIL, -1);
     CommandArgument_Int_OrDefault(MIN_TO_PRINT, 0);
     EndCommandArguments;

     // Set up output stream.

     ofstream& out = ( OUTFILE != ""
          ? *( new ofstream( OUTFILE.c_str( ) ) ) : (ofstream&) cout );

     // Parse CHUNKS.

     vec<longlong> chunks;
     if ( CHUNKS != "undefString" )
     {    int status;
          ParseLongLongSet( CHUNKS, chunks, status );
          if ( status != 0 ) ABORT( "Problem with CHUNKS argument." );    }

     // Parse MAX_FREQ.

     vec<unsigned int> max_freq;
     String MFQ = MAX_FREQ;
     while(1)
     {    if ( MFQ.IsInt( ) )
          {    max_freq.push_back( MFQ.Int( ) );
               break;    }
          if ( MFQ.Contains( ":" ) )
          {    String precolon = MFQ.Before( ":" );
               if ( !precolon.IsInt( ) )
                    ABORT( "MAX_FREQ or MF is not in correct format." );
               max_freq.push_back( precolon.Int( ) );
               MFQ = MFQ.After( ":" );    }
          else ABORT( "MAX_FREQ or MF is not in correct format." );    }

     // Parse heuristics file.

     vec< vec< pair<String, String> > > heuristics;
     vec<String> heuristics_lines;
     unsigned int max_K = K, max_MAX_FREQ = Max(max_freq);
     vec<String> allowed_heuristics;
     allowed_heuristics.push_back("K");
     allowed_heuristics.push_back("MAX_OFFSET_DIFF");
     allowed_heuristics.push_back("MIN_OVERLAP");
     allowed_heuristics.push_back("MAX_FREQ");
     allowed_heuristics.push_back("WINNING_EDGE");
     allowed_heuristics.push_back("MIN_COVERAGE");
     allowed_heuristics.push_back("MIN_HITS_TO_OVERRIDE");
     allowed_heuristics.push_back("PROGRESSION_RATIO");
     allowed_heuristics.push_back("MIN_MUTMER_LENGTH");
     allowed_heuristics.push_back("KEEP_BEST");
     allowed_heuristics.push_back("SINGLETON_HITS");
     allowed_heuristics.push_back("END_STRETCH");
     allowed_heuristics.push_back("MIN_MATCHES_PERCENT");
     allowed_heuristics.push_back("MAX_ERROR_PERCENT");
     allowed_heuristics.push_back("MIN_BASES_COVERED");

     if ( HEURISTICS != "" )
     {    if ( MAX_FREQ.Contains( ":" ) )
          {    ABORT( "Compound form of MAX_FREQ option cannot be used with "
                    << "HEURISTICS option." );    }
          Ifstream( in, HEURISTICS );
          String line, arg, left, right;
          while(1)
          {    getline( in, line );
               if ( !in ) break;
               if ( line.size( ) == 0 || line[0] == '#' ) continue;
               vec< pair<String, String> > h;
               istrstream iline( line.c_str( ) );
               while(iline)
               {    iline >> arg;
                    if ( arg.empty( ) ) break;
                    if ( !arg.Contains( "=" ) )
                    {    ABORT( "Problem in HEURISTICS file: no \"=\" found in "
                              << arg << "." );    }
                    left = arg.Before( "=" ), right = arg.After( "=" );
                    if ( left == "MO" ) left = "MAX_OFFSET_DIFF";
                    if ( left == "MF" ) left = "MAX_FREQ";
                    if ( left == "WE" ) left = "WINNING_EDGE";
                    if ( left == "MC" ) left = "MIN_COVERAGE";
                    if ( left == "PR" ) left = "PROGRESSION_RATIO";
                    if ( left == "MM" ) left = "MIN_MUTMER_LENGTH";
                    if ( left == "KB" ) left = "KEEP_BEST";
                    if ( left == "SH" ) left = "SINGLETON_HITS";
                    #define ILLHEUR(ARG)                                 \
                    {    ABORT( "Illegal value for " << ARG << " in "    \
                              << "HEURISTICS file." );    }
                    if ( left == "MAX_OFFSET_DIFF" && !right.IsInt( ) )
                         ILLHEUR( "MAX_OFFSET_DIFF or MO" );
                    if ( left == "MIN_OVERLAP" && !right.IsInt( ) )
                         ILLHEUR( "MIN_OVERLAP" );
                    if ( left == "END_STRETCH" && !right.IsInt( ) )
                         ILLHEUR( "END_STRETCH" );
                    if ( left == "MIN_MATCHES_PERCENT" && !right.IsInt( ) )
                         ILLHEUR( "MIN_MATCHES_PERCENT" );
                    if ( left == "MAX_ERROR_PERCENT" && !right.IsInt( ) )
                         ILLHEUR( "MAX_ERROR_PERCENT" );
                    if ( left == "MIN_BASES_COVERED" && !right.IsInt( ) )
                         ILLHEUR( "MIN_BASES_COVERED" );
                    if ( left == "K" && !right.IsInt( ) ) ILLHEUR( "K" );
                    if ( left == "K" )
                         max_K = Max( max_K, (unsigned int) right.Int( ) );
                    if ( left == "KEEP_BEST" && !right.IsInt( ) )
                         ILLHEUR( "KEEP_BEST or KB" );
                    if ( left == "MAX_FREQ" && !right.IsInt( ) )
                         ILLHEUR( "MAX_FREQ or MF" );
                    if ( left == "MIN_HITS_TO_OVERRIDE" && !right.IsInt( ) )
                         ILLHEUR( "MIN_HITS_TO_OVERRIDE" );
                    if ( left == "MAX_FREQ" )
                         max_MAX_FREQ
                              = Max( max_MAX_FREQ, (unsigned int) right.Int( ) );
                    if ( left == "MIN_MUTMER_LENGTH" && !right.IsInt( ) )
                         ILLHEUR( "MIN_MUTMER_LENGTH or MM" );
                    if ( left == "SINGLETON_HITS"
                         && right != "False" && right != "True" )
                    {    ILLHEUR( "SINGLETON_HITS or SH" );    }
                    double answer;
                    char c;
                    if ( left == "PROGRESSION_RATIO" )
                    {    if ( sscanf( right.c_str( ), "%lf%c", &answer, &c ) != 1 )
                              ILLHEUR( "PROGRESSION_RATIO or PR" );    }
                    if ( left == "WINNING_EDGE" )
                    {    if ( sscanf( right.c_str( ), "%lf%c", &answer, &c ) != 1 )
                              ILLHEUR( "WINNING_EDGE or WE" );    }
                    if ( left == "MIN_COVERAGE" )
                    {    if ( sscanf( right.c_str( ), "%lf%c", &answer, &c ) != 1 )
                              ILLHEUR( "MIN_COVERAGE or MC" );    }
                    if ( !Member( allowed_heuristics, left ) )
                         ABORT( "HEURISTICS file contains unknown heuristic "
                              << "parameter: " << left << "." );
                    h.push_back( make_pair( left, right ) );    }
               Sort(h);
               for ( int i = 1; i < h.isize( ); i++ )
               {    if ( h[i].first == h[i-1].first )
                         ABORT( "Same heuristic parameter appears twice in "
                              << "\"" << line << "\"" << "." );    }
               if ( h.size( ) > 0 )
               {    heuristics_lines.push_back(line);
                    heuristics.push_back(h);    }    }    }

     // Do some sanity checks on command line arguments.

     if ( TARGET_NAMING != "numeric" && TARGET_NAMING != "from_file"
          && TARGET_NAMING != "from_record"
          && TARGET_NAMING != "from_record_short"
          && TARGET_NAMING != "from_record_quoted" )
          ABORT( "TARGET_NAMING must either be \"numeric\" or \"from_file\""
               << " or \"from_record\"" << " or \"from_record_short\""
               << " or \"from_record_quoted\"." );
     if ( QUERY_NAMING != "numeric" && QUERY_NAMING != "from_record"
          && QUERY_NAMING != "from_names_file" )
     {    ABORT( "QUERY_NAMING must either be \"numeric\" or \"from_record\" "
               << " or \"from_names_file\"." );    }
     if ( !IsDirectory(TMP_DIR) )
       ABORT( "Your TMP_DIR does not exist: " + TMP_DIR );
     if ( !IsRegularFile(SEQS) )
       ABORT( "I can't find your SEQS file: " + SEQS );
     if ( QUALS != "" && !IsRegularFile(QUALS) )
          ABORT( "I can't find your QUALS file: " + QUALS );
     if ( QUALS_G != "" && !IsRegularFile(QUALS_G) )
          ABORT( "I can't find your QUALS_G file: " + QUALS_G );
     if ( !IsRegularFile(LOOKUP_TABLE) )
          ABORT( "I can't find your LOOKUP_TABLE file: " + LOOKUP_TABLE );
     if ( SEQS.Contains( ".fastb", -1 ) )
     {    if ( QUERY_NAMING != "numeric" && QUERY_NAMING != "from_names_file" )
          {    ABORT( "If SEQS is a fastb file, QUERY_NAMING must be \"numeric\" "
                    << " or \"from_names_file\"." );    }    }
     if ( QUERY_NAMING == "from_names_file" )
     {    if ( SEQ_NAMES == "" )
          {    ABORT( "You've specified QUERY_NAMING=from_names_file but "
                    << "haven't given a value for SEQ_NAMES." );    }    }
     if ( PRINT_NQS || MAX_NQS_PERCENT != undefined || MAX_QUAL_SCORE != undefined )
     {    if ( QUALS == "" )
               ABORT( "If you use PRINT_NQS or MAX_NQS_PERCENT or MAX_QUAL_SCORE, "
                    << "you have to define QUALS." );    }
     if ( SEQS_TARGETS_TO_PROCESS != "undefString" && SEQS_TO_PROCESS != "undefString" )
     {    ABORT( "Only one of SEQS_TO_PROCESS and SEQS_TARGETS_TO_PROCESS "
               << "may be specified." );    }
     if ( SEQS_TARGETS_TO_PROCESS != "undefString"
          && TARGETS_TO_PROCESS != "undefString" )
     {    ABORT( "Only one of TARGETS_TO_PROCESS and SEQS_TARGETS_TO_PROCESS "
               << "may be specified." );    }
     if ( QUALS == "" && PRINT_QUAL_SCORE )
     {    ABORT( "If PRINT_QUAL_SCORE=True, you must specify QUALS." );    }
     if ( QUALS == "" && PRINT_PREDICTED_ERRORS )
     {    ABORT( "If PRINT_PREDICTED_ERRORS=True, you must specify QUALS." );    }
     if ( QUALS_G == "" && PRINT_QUAL_SCORE )
     {    ABORT( "If PRINT_QUAL_SCORE=True, you must specify QUALS_G." );    }

     vec<fastavector> target_fasta;
     if ( TARGET_FASTA != "" ) LoadFromFastaFile( TARGET_FASTA, target_fasta );

     // Process SEQS_TO_PROCESS, SEQS_TO_EXCLUDE, TARGETS_TO_PROCESS, and
     // SEQS_TARGETS_TO_PROCESS.

     vec<longlong> seqs_to_process, to_exclude, targets_to_process;
     if ( SEQS_TO_PROCESS != "undefString" )
     {    if ( SEQS_TO_PROCESS.Contains( "random:", 0 ) )
          {    int n = SEQS_TO_PROCESS.After( "random:" ).Int( );
               ForceAssert( SEQS.Contains( ".fastb", -1 ) || SEQS_IS_FASTB );
               int N = MastervecFileObjectCount(SEQS);
               vec<Bool> keep( N, False );
               int nkeep = 0;
               while( nkeep < n )
               {    int x = randomx( ) % N;
                    if ( keep[x] ) continue;
                    keep[x] = True;
                    ++nkeep;    }
               for ( int i = 0; i < N; i++ )
                    if ( keep[i] ) seqs_to_process.push_back(i);    }
          else
          {    int status;
               ParseLongLongSet( SEQS_TO_PROCESS, seqs_to_process, status );
               if ( status != 0 )
                    ABORT( "Problem with SEQS_TO_PROCESS argument." );    }    }
     if ( SEQS_TO_EXCLUDE != "undefString" )
     {    int status;
          ParseLongLongSet( SEQS_TO_EXCLUDE, to_exclude, status );
          if ( status != 0 ) ABORT( "Problem with SEQS_TO_EXCLUDE argument." );    }
     if ( TARGETS_TO_PROCESS != "undefString" )
     {    int status;
          ParseLongLongSet( TARGETS_TO_PROCESS, targets_to_process, status );
          if ( status != 0 )
               ABORT( "Problem with TARGETS_TO_PROCESS argument." );    }
     vec< vec<longlong> > targets_for_seq;
     if ( SEQS_TARGETS_TO_PROCESS != "undefString" )
     {    vec< pair<longlong,longlong> > st;
          if ( !IsRegularFile(SEQS_TARGETS_TO_PROCESS) )
               ABORT( "Can't find file specified for SEQS_TARGETS_TO_PROCESS." );
          Ifstream( in, SEQS_TARGETS_TO_PROCESS );
          while(1)
          {    int query, target;
               in >> query;
               if ( !in ) break;
               in >> target;
               if ( !in ) ABORT( "Problem with SEQS_TARGETS_TO_PROCESS file." );
               st.push_back( make_pair( query, target ) );    }
          UniqueSort(st);
          for ( int i = 0; i < st.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < st.isize( ); j++ )
                    if ( st[j].first != st[i].first ) break;
               seqs_to_process.push_back( st[i].first );
               static vec<longlong> t;
               t.clear( );
               for ( int k = i; k < j; k++ )
               {    t.push_back( st[k].second );
                    targets_to_process.push_back( st[k].second );    }
               targets_for_seq.push_back(t);
               i = j - 1;    }
          UniqueSort(targets_to_process);    }

     // Define correspondence between query sequences in memory and their
     // numerical ids.

     // TODO: potentially dangerous truncation of index by to_query_id
     map_longlong_longlong to_query_id(IDENTITY);
     if ( SEQS_TO_PROCESS != "undefString"
          || SEQS_TARGETS_TO_PROCESS != "undefString" )
     {    to_query_id.Def( APPLY_EXTERNAL_VEC, seqs_to_process );    }

     // Open output file if needed.

     ofstream* p_outStrm = ( OUTPUT.empty() ? 0 : new ofstream( OUTPUT.c_str() ) );

     // Read in sequences.

     vecbasevector seq;
     vec<fastavector> seq_fasta;
     bool fastb = false;
     if ( SEQS.Contains( ".fastb", -1 ) || SEQS_IS_FASTB ) {
       fastb = true;
       if ( !IsRegularFile(SEQS) )
	 ABORT( SEQS << " does not exist." );
       if ( !IsGoodFeudalFile(SEQS) ) {
	 ABORT( SEQS << " is not in fastb format." );
       }
     }

     //Simple case: all sequences.
     if ( "undefString" ==  SEQS_TO_PROCESS &&
	  "undefString" == SEQS_TARGETS_TO_PROCESS ) {
       if (fastb) seq.ReadAll(SEQS);
       else 
       {    FetchReads( seq, 0, SEQS );
            LoadFromFastaFile( SEQS, seq_fasta );    }
     }
     else { //Read in only selected sequences

       longlong maxSeqs=0;
       if (fastb) maxSeqs = MastervecFileObjectCount(SEQS);
       else maxSeqs = LineOfOutput("grep -c \">\" " + SEQS).Int();

       //Shorten to maxSeqs if possible, abort otherwise.
       if ( seqs_to_process.nonempty( )
	    && seqs_to_process.back( ) >= maxSeqs ) {
	 if (seqs_to_process.front() >= maxSeqs) {
	   ABORT( "You've requested query sequence number "
		  << seqs_to_process.back( )
		  << ".  There aren't that many\nsequences." );
	 }
	 else { //use all possible sequences: remove bad ones.
	   seqs_to_process.erase
	     (remove_if(seqs_to_process.begin(), seqs_to_process.end(),
			bind1st(less_equal<int>(), maxSeqs) ),
	      seqs_to_process.end());
	 }
       }

       if (fastb) seq.Read( SEQS, seqs_to_process );
       else 
       {    // Note truncation here....
            vec<int> seqs_to_process_int;
            for ( size_t i = 0; i < seqs_to_process.size( ); i++ )
                 seqs_to_process_int.push_back( seqs_to_process[i] );
            FetchReads( seq, 0, SEQS, 0, 0, out, &seqs_to_process_int );
            vec<fastavector> seq_fasta_all;
            LoadFromFastaFile( SEQS, seq_fasta_all );
            for ( int i = 0; i < seqs_to_process.isize( ); i++ )
                 seq_fasta.push_back( seq_fasta_all[ seqs_to_process[i] ] );    }
     }

     // Define sequence names.

     vec<String> seq_names;
     seq_names.reserve( seq.size( ) );
     if ( QUERY_NAMING == "numeric" )
     {    for ( size_t i = 0; i < seq.size( ); i++ )
               seq_names.push_back( ToString( to_query_id(i) ) );    }
     else if ( QUERY_NAMING == "from_record" )
     {
       if ( SEQS.Contains( ".fastb", -1 ) )
	 ABORT( "QUERY_NAMING=from_record not possible for a fastb file." );

       fast_ifstream nin(SEQS);
       String line;
       longlong count = 0;
       while(1)
       {    getline( nin, line );
            if ( nin.fail( ) ) break;
            if ( line.Contains( ">", 0 ) )
            {    if ( "undefString" != SEQS_TO_PROCESS
                      && !BinMember( seqs_to_process, count++ ) )
                 {    continue;    }
                 if ( line.Contains( " " ) )
                      seq_names.push_back( line.After( ">" ).Before( " " ) );
                 else seq_names.push_back( line.After( ">" ) );    }    }

       /*
       FastaNameParser *name_parser;
       if ( nameParser == "first" ) name_parser = new FirstWordParser;
       else name_parser = new LastWordParser;
       FastaSequenceFilestream filestream_to_read( SEQS, name_parser );
       vecString all_names;
       filestream_to_read.getOnlyNames( all_names );
       for ( int i = 0; i < all_names.size(); ++i )
         seq_names.push_back( all_names[i] );
       delete name_parser;
       */
     }
     else if ( QUERY_NAMING == "from_names_file" )
     {    fast_ifstream in(SEQ_NAMES);
          String line;
          int count = 0;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( ( SEQS_TO_PROCESS == "undefString"
                         && SEQS_TARGETS_TO_PROCESS == "undefString" )
                    || BinMember(seqs_to_process, count) )
                    seq_names.push_back(line);
               ++count;    }    }
     if ( seq.size( ) != seq_names.size( ) )
     {    ABORT( "The number of query sequences (" << seq.size( ) << ") does not "
               << "agree with the number of query\nsequence names "
               << "(" << seq_names.size( ) << ")." );    }
     for ( size_t i = 0; i < seq.size( ); i++ )
          if ( seq[i].size( ) > max_query_size )
          {    ABORT( "One of your query sequences has length " << seq[i].size( )
                    << ".  The maximum allowed length is 200 Mb." );    }

     // Set up sequence mask.  The idea is that if we find a partial alignment,
     // and no full alignment, then on subsequence major passes, we look for
     // alignments which are disjoint from the partial alignment.

     vec< vec< pair<int,int> > > seq_mask( seq.size( ) );

     // Read in quality scores for sequences, if given.

     vecqualvector qual;
     if ( QUALS.Contains( ".qualb", -1 ) || QUALS_IS_QUALB )
     {    if ( !IsGoodFeudalFile(QUALS) )
          {    ABORT( QUALS << " is not in qualb format." );    }
          if ( SEQS_TO_PROCESS != "undefString"
                    || SEQS_TARGETS_TO_PROCESS != "undefString" )
          {    if ( seqs_to_process.nonempty( )
                    && static_cast<size_t>(seqs_to_process.back()) >= MastervecFileObjectCount(QUALS) )
               {    ABORT( "You've requested query sequence number "
                         << seqs_to_process.back( )
                         << ".  There aren't that many\nqual sequences." );    }
               qual.Read( QUALS, seqs_to_process );    }
          else qual.ReadAll(QUALS);    }
     else if ( QUALS != "" )
     {    if ( SEQS_TO_PROCESS != "undefString"
                    || SEQS_TARGETS_TO_PROCESS != "undefString" )
          {    vec<int> seqs_to_process_int;
               for ( size_t i = 0; i < seqs_to_process.size( ); i++ )
                    seqs_to_process_int.push_back( seqs_to_process[i] );
               ReadFastaQuals( QUALS, qual, &seqs_to_process_int );    }
          else ReadFastaQuals( QUALS, qual );    }

     // Truncate reads.

     if ( TRUNCATE_TO >= 0 )
     {    for ( size_t i = 0; i < seq.size( ); i++ )
          {    if ( seq[i].isize( ) > TRUNCATE_TO )
                    seq[i].resize(TRUNCATE_TO);
               if ( qual.size( ) > 0 && (int) qual[i].size( ) > TRUNCATE_TO )
                    qual[i].resize(TRUNCATE_TO);    }    }
     if ( TRUNCATE_TO_TAIL >= 0 )
     {    for ( size_t i = 0; i < seq.size( ); i++ )
          {    if ( seq[i].isize( ) > TRUNCATE_TO_TAIL )
               {    basevector b = seq[i];
                    seq[i].resize(TRUNCATE_TO_TAIL);    
                    for ( int j = 0; j < TRUNCATE_TO_TAIL; j++ )
                         seq[i].Set( j,  
                              b[ b.size( ) - TRUNCATE_TO_TAIL + j ] );    }    }    }

     // Read in quality scores for genome, if given.

     vecqualvector qual_g;
     if ( QUALS_G.Contains( ".qualb", -1 ) )
     {    if ( !IsGoodFeudalFile(QUALS_G) )
          {    ABORT( QUALS_G << " is not in qualb format." );    }
          qual_g.ReadAll(QUALS_G);    }
     else if ( QUALS_G != "" ) ReadFastaQuals( QUALS_G, qual_g );

     // Do sanity check on consistency of seqs and quals.

     if ( qual.size( ) )
     {    if ( seq.size( ) != qual.size( ) )
          {    ABORT( "You've supplied " << seq.size( ) << " query sequences "
                    << "but " << qual.size( ) << " quality score sequences." );    }
          for ( size_t i = 0; i < seq.size( ); i++ )
          {    if ( seq[i].size( ) != qual[i].size( ) )
               {    ABORT( "Query sequence " << i << " has length "
                         << seq[i].size( ) << ", but you've supplied "
                         << qual[i].size( ) << " quality scores." );    }    }    }

     // Read header information from lookup table.
     lookup_table look(LOOKUP_TABLE);

     unsigned int Ktab = look.K( );
     if ( Ktab < max_K )
     {    ABORT( "You've specified a K value of " << max_K << " on the command "
               << "line or in the HEURISTICS file, but the K value defined by the "
               << "lookup table is smaller: " << Ktab << "." );    }

     // Set up hits and hitsp vectors.

     vec<unsigned int>* hits_ptr = 0;
     vec< pair<unsigned int, unsigned int> >* hitsp_ptr = 0;

     // Set up for global statistics.

     vec<int> small_insertions(1000), large_insertions;
     vec<int> small_deletions(1000), large_deletions;
     longlong mismatches = 0, total_aligned = 0;

     // Define target names.

     vec<String> target_names( look.NContigs( ) );
     for ( unsigned int i = 0; i < look.NContigs( ); i++ )
     {    static String targ;
          if ( TARGET_NAMING == "numeric" ) target_names[i] = ToString(i);
          else if ( TARGET_NAMING == "from_file" )
               target_names[i] = look.ContigNameBasic(i);
          else if ( TARGET_NAMING == "from_record" )
               target_names[i] = look.ContigName(i);
          else if ( TARGET_NAMING == "from_record_short" )
          {    targ = look.ContigName(i);
               int j;
               for ( j = 0; j < (int) targ.size( ); j++ )
                    if ( isspace( targ[j] ) ) break;
               target_names[i] = targ.substr( 0, j );    }
          else
          {    targ = look.ContigName(i);
               targ.GlobalReplaceBy( "\"", "" );
               target_names[i] = "\"" + targ + "\"";    }    }

     // Go through the lookup table chunk by chunk.

     vec<Bool> placed( seq.size( ), False ), accepted( seq.size( ), False );
     vec<Bool> placed_full_length( seq.size( ), False );
     vec<Bool> uniquely_placed( seq.size( ), False );
     vec<Bool> progress( seq.size( ), False );
     int major_pass_count = max_freq.size( );
     if ( heuristics.size( ) > 0 ) major_pass_count = heuristics.size( );
     for ( int major_pass = 0; major_pass < major_pass_count; major_pass++ )
     {
          // Announce pass.

          if ( !QUIET )
          {    out << "\n" << Date( ) << ": STARTING MAJOR PASS " << major_pass + 1
                    << endl;    }
          Bool same_heuristics = ( heuristics.nonempty( )
               && major_pass >= 1
               && heuristics[major_pass-1] == heuristics[major_pass] );
          if ( heuristics.size( ) > 0 )
          {    out << "(from HEURISTICS:";
               istrstream ih( heuristics_lines[major_pass].c_str( ) );
               String opt;
               while(ih)
               {    ih >> opt;
                    if ( opt.empty( ) ) break;
                    out << " " << opt;    }
               out << ")" << endl;    }

          // Define a macro "LOCALIZE" which allows a variable to be redefined
          // in the current scope and set equal to the same variable, defined in
          // the parent scope.

          #define LOCALIZE( TYPE, NAME )            \
               TYPE NAME ## _localize_temp = NAME;  \
               TYPE NAME = NAME ## _localize_temp;

          // Localize and handle K.

          LOCALIZE( unsigned int, K );
          if ( heuristics.size( ) > 0 )
          {    vec< pair<String, String> > heur = heuristics[major_pass];
               for ( int i = 0; i < heur.isize( ); i++ )
               {    if ( heur[i].first == "K" )
                         K = heur[i].second.Int( );    }    }
          unsigned int Kdiffbits = 2 * (Ktab - K);
          unsigned int four_to_Kdiff = 1;
          for ( unsigned int i = K; i < Ktab; i++ )
               four_to_Kdiff *= 4;

          // Localize WINNING_EDGE.

          LOCALIZE( double, WINNING_EDGE );
          if ( heuristics.size( ) > 0 )
          {    vec< pair<String, String> > heur = heuristics[major_pass];
               for ( int i = 0; i < heur.isize( ); i++ )
               {    if ( heur[i].first == "WINNING_EDGE" )
                    {    char c;
                         sscanf( heur[i].second.c_str( ), "%lf%c",
                              &WINNING_EDGE, &c );    }    }    }

          // Localize PROGRESSION_RATIO.

          LOCALIZE( double, PROGRESSION_RATIO );
          if ( heuristics.size( ) > 0 )
          {    vec< pair<String, String> > heur = heuristics[major_pass];
               for ( int i = 0; i < heur.isize( ); i++ )
               {    if ( heur[i].first == "PROGRESSION_RATIO" )
                    {    char c;
                         sscanf( heur[i].second.c_str( ), "%lf%c",
                              &PROGRESSION_RATIO, &c );    }    }    }

          // Localize MIN_COVERAGE.

          LOCALIZE( double, MIN_COVERAGE );
          if ( heuristics.size( ) > 0 )
          {    vec< pair<String, String> > heur = heuristics[major_pass];
               for ( int i = 0; i < heur.isize( ); i++ )
               {    if ( heur[i].first == "MIN_COVERAGE" )
                    {    char c;
                         sscanf( heur[i].second.c_str( ), "%lf%c",
                              &MIN_COVERAGE, &c );    }    }    }

          // Localize integer parameters.

          #define LOCALIZE_INT_PARAM(PARAM)                               \
               LOCALIZE( unsigned int, PARAM );                           \
               if ( heuristics.nonempty( ) )                              \
               {    vec< pair<String, String> >& heur                     \
                         = heuristics[major_pass];                        \
                    for ( int i = 0; i < heur.isize( ); i++ )             \
                    {    if ( heur[i].first == #PARAM )                   \
                              PARAM = heur[i].second.Int( );    }    }
          LOCALIZE_INT_PARAM(MIN_HITS_TO_OVERRIDE);
          LOCALIZE_INT_PARAM(MAX_OFFSET_DIFF);
          LOCALIZE_INT_PARAM(MIN_OVERLAP);
          LOCALIZE_INT_PARAM(MIN_MUTMER_LENGTH);
          LOCALIZE_INT_PARAM(KEEP_BEST);
          LOCALIZE_INT_PARAM(END_STRETCH);
	  LOCALIZE_INT_PARAM(MIN_MATCHES_PERCENT);
	  LOCALIZE_INT_PARAM(MAX_ERROR_PERCENT);
	  LOCALIZE_INT_PARAM(MIN_BASES_COVERED);

          // Localize MAX_FREQ.

          unsigned int mfreq;
          if ( heuristics.empty( ) ) mfreq = max_freq[major_pass];
          else
          {    mfreq = max_freq[0];
               vec< pair<String, String> > heur = heuristics[major_pass];
               for ( int i = 0; i < heur.isize( ); i++ )
               {    if ( heur[i].first == "MAX_FREQ" )
                         mfreq = heur[i].second.Int( );    }    }

          // Localize and handle SINGLETON_HITS.

          LOCALIZE( Bool, SINGLETON_HITS );
          if ( heuristics.size( ) > 0 )
          {    vec< pair<String, String> > heur = heuristics[major_pass];
               for ( int i = 0; i < heur.isize( ); i++ )
               {    if ( heur[i].first == "SINGLETON_HITS" )
                    {    if ( heur[i].second == "True" ) SINGLETON_HITS = True;
                         else SINGLETON_HITS = False;    }    }    }
          if (SINGLETON_HITS)
          {    if ( hitsp_ptr != 0 )
               {    delete hitsp_ptr;
                    hitsp_ptr = 0;    }
               if ( hits_ptr == 0 )
               {    hits_ptr = new vec<unsigned int>;
                    hits_ptr->reserve( 5 * max_MAX_FREQ );    }    }
          else
          {    if ( hits_ptr != 0 )
               {    delete hits_ptr;
                    hits_ptr = 0;    }
               if ( hitsp_ptr == 0 )
               {    hitsp_ptr = new vec< pair<unsigned int, unsigned int> >;
                    hitsp_ptr->reserve( 5 * max_MAX_FREQ );    }    }

          // Set up temporary files for alignments.

          String real_tmp = RealPath(TMP_DIR);
          temp_file qualifiers_file( real_tmp + "/QueryLookupTable_tmp1_XXXXXXX" );
          int fdq = OpenForWrite(qualifiers_file);
          temp_file qualifiers_summary_file(
               real_tmp + "/QueryLookupTable_tmp2_XXXXXXX" );
          int fdqs = OpenForWrite(qualifiers_summary_file);
          fchmod( fdq, 0666 ), fchmod( fdqs, 0666 );
          off_t disk_pos = 0;

          int n_all_qualifiers = 0;
          for ( unsigned int i = 0; i < look.NChunks( ); i++ )
          {
               // Read in a chunk.

               Bool ignore_chunk = ( CHUNKS != "undefString"
                    && !BinMember( chunks, (int) i+1 ) );
               if ( targets_to_process.nonempty( ) )
               {    unsigned int chunk_start = look.StartBaseInChunk(i);
                    unsigned int chunk_stop = look.StopBaseInChunk(i);
                    Bool need_chunk = False;
                    for ( int j = 0; j < targets_to_process.isize( ); j++ )
                    {    int t = targets_to_process[j];
                         unsigned int contig_start = look.ContigStart(t);
                         unsigned int contig_stop = look.ContigStop(t);
                         unsigned int over = IntervalOverlap( contig_start,
                              contig_stop, chunk_start, chunk_stop );
                         if ( over > 0 ) need_chunk = True;    }
                    if ( !need_chunk ) continue;    }
               if ( ANNOUNCE_ITERATIONS && !ignore_chunk )
               {    out << "\n" << Date( ) << ": processing chunk " << i+1 << " of "
                         << look.NChunks( ) << "\n";
                    if (PRINT_MEMORY_USAGE) PrintMemUsage(out);
                    out << endl;    }
               if ( ignore_chunk && (int) i+1 > Max(chunks) ) break;
               if ( ignore_chunk ) continue;
               look.ReadChunk(i);

               // Go through the sequences.

               for ( size_t j = 0; j < seq.size( ); j++ )
               {    if ( accepted[j] ) continue;
                    if ( same_heuristics && !progress[j] ) continue;
                    int query_id = to_query_id(j);
                    if ( BinMember( to_exclude, query_id ) ) continue;

                    // If SEQS_TARGETS_TO_PROCESS has been specified, check to see
                    // if any of the targets for this sequence are in this chunk.

                    if ( SEQS_TARGETS_TO_PROCESS != "undefString" )
                    {    Bool target_in_chunk = False;
                         for ( int l = 0; l < targets_for_seq[j].isize( ); l++ )
                         {    int t = targets_for_seq[j][l];
                              if ( IntervalOverlap( look.ContigStart(t),
                                   look.ContigStop(t), look.StartBaseInChunk(i),
                                   look.StopBaseInChunk(i) ) > 0 )
                              {    target_in_chunk = True;
                                   break;    }    }
                         if ( !target_in_chunk ) continue;    }

                    static basevector s;
                    s = seq[j];
                    if ( s.size( ) < K ) continue;
                    if ( s.size( ) < MIN_QUERY_LENGTH ) continue;
                    static qualvector q;
                    if ( qual.size( ) ) q = qual[j];

                    static vec< pair<int, int> > ints;
                    ints.clear( );

                    if ( seq_mask[j].size( ) == 0 )
                         ints.push_back( make_pair( 0, s.size( ) ) );
                    else
                    {    static vec<ho_interval> covered;
                         covered.clear( );
                         for ( int u = 0; u < seq_mask[j].isize( ); u++ )
                              covered.push_back( ho_interval(
                                   seq_mask[j][u].first, seq_mask[j][u].second ) );
                         vec< pair<ho_interval, int> > condensed;
                         CondenseIntervals( s.size( ), covered, condensed );
                         for ( int u = 0; u < condensed.isize( ); u++ )
                         {    if ( condensed[u].second != 0 ) continue;
                              const ho_interval& h = condensed[u].first;
                              if ( h.Length( ) < (int) K ) continue;
                              ints.push_back(
                                   make_pair( h.Start( ), h.Stop( ) ) );    }    }
                    for ( int u = 0; u < ints.isize( ); u++ )
                    {    const vec<longlong>& targets =
                              ( SEQS_TARGETS_TO_PROCESS == "undefString"
                                   ? targets_to_process : targets_for_seq[j] );
                         int npasses = ( FW_ONLY ? 1 : 2 );
                         Bool firstpass = ( RC_ONLY ? 2 : 1 );
                         ProcessQuerySequence( npasses, firstpass, look, j, s,
                              ints[u].first, ints[u].second, s.size( ), q,
                              n_all_qualifiers, disk_pos, Kdiffbits, four_to_Kdiff,
                              mfreq, to_query_id, K, MIN_COVERAGE,
                              MIN_HITS_TO_OVERRIDE, WINNING_EDGE, MAX_OFFSET_DIFF,
                              MIN_OVERLAP, END_STRETCH, PROGRESSION_RATIO,
                              MIN_MUTMER_LENGTH, SINGLETON_HITS, KEEP_BEST,
                              SHOW_PREMUTMERS, SHOW_MUTMERS, hits_ptr, hitsp_ptr,
                              fdq, fdqs, IMPERFECT_EXTENSION,
                              targets, SYNC_ALIGNS_TO_TACG, SMITH_WAT,
                              BW_ADD, SW_MISMATCH_PENALTY, SW_GAP_PENALTY,
                              SW_GAP_VERBOSE, out );    }    }    }

          // Free up some space.

          if ( major_pass == major_pass_count - 1 )
          {    if ( hits_ptr != 0 ) delete hits_ptr;
               if ( hitsp_ptr != 0 ) delete hitsp_ptr;
               look.DestroyStuff( );    }

          // Summarize qualifiers.

          progress.resize_and_set( seq.size( ), False );
          close(fdq);
          fdq = OpenForRead(qualifiers_file);
          if ( fdq < 0 )
          {    ABORT( "Unable to reopen alignment scratch file.  Could it have been "
                    << "deleted somehow?" );    }
          close(fdqs);
          vec< pair<longlong, off_t> > all_qualifiers_summary( n_all_qualifiers );
          if ( n_all_qualifiers > 0 )
          {    FileReader fr(qualifiers_summary_file.c_str());
               fr.read( &all_qualifiers_summary[0],
                       (sizeof(longlong) + sizeof(off_t))*n_all_qualifiers );  }
          Remove(qualifiers_summary_file);
          Sort(all_qualifiers_summary);

          if ( !QUIET )
          {    out << "\n" << Date( ) << ": RESULTS OF PASS " << major_pass + 1
                    << endl;    }
          if (PRINT_MEMORY_USAGE) PrintMemUsage(out);
          for ( int i = 0; i < all_qualifiers_summary.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < all_qualifiers_summary.isize( ); j++ )
                    if ( all_qualifiers_summary[j].first
                         != all_qualifiers_summary[i].first ) break;

               int id = all_qualifiers_summary[i].first;

               // Read in the alignments.

               static vec<look_align> these_qualifiers;
               int q_count = j - i;
               if ( q_count > these_qualifiers.isize( ) )
                    these_qualifiers.resize(q_count);
               for ( int x = 0; x < q_count; x++ )
               {    off_t disk_loc = all_qualifiers_summary[ i + x ].second;
                    lseek( fdq, disk_loc, SEEK_SET );
                    these_qualifiers[x].BinaryRead(fdq);    }
               sort( these_qualifiers.begin( ),
                    these_qualifiers.begin( ) + q_count );

               // Compute minimum number of hits required.  If all alignments have
               // less hits than the number required by MIN_COVERAGE, discard them
               // all.

               int best_hit_count = 0;
               for ( int m = i; m < j; m++ )
                    best_hit_count
                         = Max( best_hit_count, these_qualifiers[m-i].nhits );

               int holds = seq[id].size( ) - K + 1;
               if ( seq_mask[id].size( ) > 0 )
               {    holds = 0;
                    static vec<ho_interval> covered;
                    covered.clear( );
                    for ( int u = 0; u < seq_mask[id].isize( ); u++ )
                         covered.push_back( ho_interval(
                              seq_mask[id][u].first, seq_mask[id][u].second ) );
                    vec< pair<ho_interval, int> > condensed;
                    CondenseIntervals( seq[id].size( ), covered, condensed );
                    for ( int u = 0; u < condensed.isize( ); u++ )
                    {    if ( condensed[u].second != 0 ) continue;
                         const ho_interval& h = condensed[u].first;
                         holds = Max( holds, h.Length( ) - (int) K + 1 );    }    }

               int hits_floor = int( floor( MIN_COVERAGE * ((double) holds) ) );
               hits_floor = Min( hits_floor, (int) MIN_HITS_TO_OVERRIDE );
               int winning_edge
                    = int( floor( WINNING_EDGE * sqrtf(float(best_hit_count)) ) );
               int min_hits = Max( 0, best_hit_count - winning_edge );
               if ( best_hit_count < hits_floor )
               {    i = j - 1;
                    continue;    }

               // If an alignment is dominated by another alignment, mark it for
               // deletion.  If there is a full alignment, delete any partial
               // alignment whose error rate is 50% greater than the full alignment.
               // If there are two highly overlapping alignments, try to pick one.

               float full_error_rate = 1000000.0;
               for ( int k = i; k < j; k++ )
               {    const look_align& l = these_qualifiers[ k - i ];
                    if ( l.FullLength( ) )
                         full_error_rate
                              = Min( full_error_rate, l.ErrorRate( ) );    }
               static vec<Bool> dominated;
               dominated.resize_and_set( j - i, False );
               if (REMOVE_DOMINATED)
               {    for ( int k1 = i+1; k1 < j; k1++ )
                    {    const look_align& l1 = these_qualifiers[ k1 - i ];
                         if ( !l1.FullLength( )
                              && l1.ErrorRate( ) >= 1.5 * full_error_rate )
                         {    dominated[k1-i] = True;    }
                         for ( int k2 = i; k2 < j; k2++ )
                         {    const look_align& l2 = these_qualifiers[ k2 - i ];
                              if ( l1.target_id != l2.target_id ) continue;
                              if ( l1.a.pos2( ) > l2.a.pos2( )
                                   && l1.a.Pos2( ) <= l2.a.Pos2( ) )
                              {    dominated[k1-i] = True;    }
                              if ( l1.a.pos2( ) >= l2.a.pos2( )
                                   && l1.a.Pos2( ) < l2.a.Pos2( ) )
                              {    dominated[k1-i] = True;    }
                              if ( k1 > k2 && l1.a.pos1( ) == l2.a.pos1( )
                                   && l1.a.Pos1( ) == l2.a.Pos1( )
                                   && l1.a.pos2( ) == l2.a.pos2( )
                                   && l1.a.Pos2( ) == l2.a.Pos2( )
                                   && l1.mutations >= l2.mutations )
                              {    dominated[k1-i] = True;    }
                              if ( l1.a.pos1( ) == l2.a.pos1( )
                                   && l1.a.Pos1( ) == l2.a.Pos1( )
                                   && float( IntervalOverlap( l1.a.pos2( ),
                                        l1.a.Pos2( ), l2.a.pos2( ), l2.a.Pos2( ) ) )
                                        >= 0.8 * float(
                                        l1.a.Pos2( ) - l1.a.pos2( ) ) )
                              {    if ( l1.mutations > l2.mutations )
                                        dominated[k1-i] = True;    }    }    }    }

               // The alignments of query id are i..j.  If the first has full length,
               // print only full length alignments having <= 2 times as many errors
               // as the first alignment or <= 8 more errors than the first
               // alignment.  Note that since the first alignment is not necessarily
               // the best (since alignments were sorted only by mutations), this
               // doesn't make total sense.  Also, filter based on number of hits.
               //
               // Also, print no more than KEEP_BEST alignments.

               // First count the number of alignments which we would consider
               // printing.

               unsigned int to_print = 0;
               for ( int m = i; m < j; m++ )
               {    const look_align& q = these_qualifiers[m-i];

                    if (FILTER)
                    {    if ( dominated[ m - i ] ) continue;
                         if ( m > i && q.nhits < min_hits ) continue;
                         if ( these_qualifiers[0].FullLength( )
                              && q.Errors( ) > 2 * these_qualifiers[0].Errors( )
                              && q.Errors( )
                                   > FILTER_ADD + these_qualifiers[0].Errors( ) )
                         {    continue;    }    }
                    ++to_print;    }
               if ( (int) to_print < MIN_TO_PRINT )
               {    i = j - 1;
                    continue;    }

               // Process the alignments.

               int print_count = 0;
               for ( int m = i; m < j; m++ )
               {    const look_align& q = these_qualifiers[m-i];

                    // See if the alignment is acceptable.

                    if (FILTER)
                    {    if ( dominated[ m - i ] ) continue;
                         if ( m > i && q.nhits < min_hits ) continue;
                         if ( these_qualifiers[0].FullLength( )
                              && q.Errors( ) > FILTER_MULT * these_qualifiers[0].Errors( )
                              && q.Errors( )
                                   > FILTER_ADD + these_qualifiers[0].Errors( ) )
                         {    continue;    }    }

                    // Check if over MAX_PLACEMENTS.

                    if ( MAX_PLACEMENTS != undefined && to_print > MAX_PLACEMENTS )
                         continue;

                    // Fetch genome bases.

                    unsigned int start = (unsigned int) -1;
                    static basevector t;
                    if ( PARSEABLE || VISUAL || PRINT_NQS || PRINT_RMR
                         || RMR_BY_BLOCK || MAX_RMR_PERCENT != undefined
                         || MIN_MM_PERCENT != undefined
                         || MAX_NQS_PERCENT != undefined
                         || MAX_QUAL_SCORE != undefined || PRINT_QUAL_SCORE )
                    {    unsigned int begin = look.ContigStart( q.target_id );
                         start = 0;
                         if ( 80 < q.a.pos2( ) ) start = q.a.pos2( ) - 80;
                         unsigned int stop = q.a.Pos2( ) + 80;
                         if ( !look.CanFetchBasesFromDisk(
                              begin + start, begin + stop ) )
                         {    start = q.a.pos2( );
                              stop = q.a.Pos2( );    }
                         t.Setsize( stop - start );
                         look.FetchBasesFromDisk( begin + start,
                              begin + stop, t );    }

                    // Check for improper alignments.

                    if (REQUIRE_PROPER)
                    {    if ( q.a.pos1( ) > 0 && q.a.pos2( ) > 0 ) continue;
                         if ( (unsigned int) q.a.Pos1( ) < q.query_length
                              && (unsigned int) q.a.Pos2( ) < q.target_length )
                         {    continue;    }    }
                    if (REQUIRE_FULL1)
                    {    if ( q.a.pos1( ) > 0 ) continue;
                         if ( (unsigned int) q.a.Pos1( ) < q.query_length )
                         {    continue;    }    }
                    if (REQUIRE_FULL2)
                    {    if ( q.a.pos2( ) > 1 ) continue;
                         if ( (unsigned int) q.a.Pos2( ) < q.target_length - 1 )
                         {    continue;    }    }
                    if (REQUIRE_POS2_ZERO)
                    {    if ( q.a.pos2( ) > 0 ) continue;    }

                    // Check for too many mismatches or indels.

                    if ( q.mutations > (int) MAX_MISMATCHES ) continue;
                    if ( q.indels > (int) MAX_INDELS ) continue;
		    if ( q.ErrorRate()*100 > (int)MAX_ERROR_PERCENT) continue;
		    // Check for too many mismatches/alignment length.
		    FlowAlignSummary fas(q);
		    if ( fas.coverScorePercent() < int(MIN_MATCHES_PERCENT))
		      continue;
		    // Check for too short an alignment.
		   if (fas.goodBases() < int(MIN_BASES_COVERED)) continue;

                    // Check for too many high quality discrepancies or too high
                    // a reciprocal mismatch rate.

                    if ( MIN_PERFECT_MATCH > 0 )
                    {    basevector s = seq[q.query_id];
                         if ( q.Rc1( ) ) s.ReverseComplement( );
                         vec<ho_interval> perfs;
                         q.a.PerfectIntervals1( s, t, perfs );
                         int min_perf = 0;
                         for ( int j = 0; j < perfs.isize( ); j++ )
                              min_perf = Max( min_perf, perfs[j].Length( ) );
                         if ( min_perf < MIN_PERFECT_MATCH ) continue;    }
                    if ( MAX_NQS_PERCENT != undefined )
                    {    int see = q.CountNqs( seq[q.query_id], qual[q.query_id],
                              t, start );
                         if ( 100.0 * float(see) / float( q.a.Pos1( ) - q.a.pos1( ) )
                              > MAX_NQS_PERCENT )
                         {    continue;    }    }
                    if ( MAX_RMR_PERCENT != undefined )
                    {    if ( q.RmrPercent( seq[q.query_id], qual[q.query_id],
                              t, start ) > Float(MAX_RMR_PERCENT) )
                         {    continue;    }    }
                    if ( MIN_MM_PERCENT != undefined )
                    {    if ( q.MMPercent( seq[q.query_id], qual[q.query_id],
                              t, start ) < Float(MIN_MM_PERCENT) )
                         {    continue;    }    }
                    if ( MAX_QUAL_SCORE != undefined )
                    {    static qualvector qnull(0);
                         const qualvector& qt =
                              ( QUALS_G != "" ? qual_g[q.target_id] : qnull );
                         if ( q.QualScore( seq[ q.query_id ],
                              qual[ q.query_id ], t, start, qt )
                              > Float(MAX_QUAL_SCORE) )
                         {    continue;    }    }

                    // Update sequence mask.

                    if ( !these_qualifiers[0].FullLength( ) )
                    {    if ( !q.rc1 )
                         {    seq_mask[ q.query_id ].push_back(
                                   make_pair( q.a.pos1( ), q.a.Pos1( ) ) );    }
                         else
                         {    seq_mask[ q.query_id ].push_back(
                                   make_pair( q.query_length - q.a.Pos1( ),
                                        q.query_length - q.a.pos1( ) ) );    }    }
                    else seq_mask[ q.query_id ].clear( );
                    progress[ q.query_id ] = True;

                    // Record global statistics.

                    if (GLOBAL_STATS)
                    {    const basevector &query = seq[ q.query_id ], &target = t;
                         total_aligned += q.a.Pos1( ) - q.a.pos1( );
                         int p1 = q.a.pos1( ), p2 = q.a.pos2( );
                         for ( int j = 0; j < q.a.Nblocks( ); j++ )
                         {    int g = q.a.Gaps(j);
                              if ( g > 0 )
                              {    p2 += g;
                                   if ( g < 1000 ) ++small_deletions[g];
                                   else large_deletions.push_back(g);    }
                              if ( g < 0 )
                              {    p1 -= g;
                                   if ( -g < 1000 ) ++small_insertions[-g];
                                   else large_insertions.push_back(-g);    }
                              for ( int x = 0; x < q.a.Lengths(j); x++ )
                              {    Bool mismatch = False;
                                   if ( !q.rc1 && query[p1] != target[ p2 - start ] )
                                   {    mismatch = True;    }
                                   if ( q.rc1 && 3 - query[q.query_length - p1 - 1]
                                        != target[ p2 - start ] )
                                   {    mismatch = True;    }
                                   if (mismatch) ++mismatches;
                                   ++p1, ++p2;    }    }    }

                    // Print.

                    Bool print_something = PARSEABLE || PARSEABLE_BRIEF
                         || READABLE_BRIEF || PRINT_RMR || PRINT_MM || RMR_BY_BLOCK
                         || PRINT_NQS || PRINT_QUAL_SCORE || VISUAL
                         || DUMP_WINDOW > 0 || PRINT_PREDICTED_ERRORS;
		    ostream *outp = OUTPUT.empty( ) ? ( ostream* ) &out : p_outStrm;
                    ostream &out = *outp;
		    if ( print_count == 0 && print_something ) out << "\n";
                    static qualvector qnull(0);
                    const qualvector& qid =
                         ( qual.size( ) ? qual[ q.query_id ] : qnull );
                    if (PARSEABLE)
                         q.PrintParseable( out, seq[ q.query_id ], qid,
                              t, start, seq_names, target_names );
                    if (PARSEABLE_BRIEF)
                         q.PrintParseableBrief( out, seq[ q.query_id ], qid,
                              t, start, seq_names, target_names );
                    if (READABLE_BRIEF)
                         q.PrintReadableBrief( out, seq[ q.query_id ], qid,
                              t, start, seq_names, target_names );
                    if (PRINT_PREDICTED_ERRORS)
                    {    double errs = 0.0;
                         for ( qvec::size_type i = 0; i < qid.size( ); i++ )
                              errs += pow( 10.0, -double(qid[i])/10.0 );
                         out << "predicted errors = " << errs << "\n";    }
                    if (PRINT_RMR)
                         q.PrintRmr( out, seq[ q.query_id ], qid, t, start );
                    if (PRINT_MM)
                         q.PrintMM( out, seq[ q.query_id ], qid, t, start );
                    if (RMR_BY_BLOCK)
                         q.PrintRmrByBlock( out, seq[ q.query_id ], qid,
                              t, start, seq_names, target_names );
                    if (PRINT_NQS)
                         q.PrintNqs( out, seq[ q.query_id ], qual[ q.query_id ],
                              t, start, seq_names, target_names );
                    if (PRINT_QUAL_SCORE)
                         q.PrintQualScore( out, seq[ q.query_id ],
                              qual[ q.query_id ], t, start, qual_g[q.target_id] );
                    if (VISUAL)
                    {    if ( QUALS_G == "" )
                         {    if ( TARGET_FASTA.size( ) == 0 )
                              {    if ( seq_fasta.size( ) == 0 )
                                   {    q.PrintVisual( out, seq[ q.query_id ], qid,
                                             t, start, VISUAL_ABBR );    }
                                   else
                                   {    q.PrintVisual( out, seq_fasta[ q.query_id ],
                                             qid, t, start, VISUAL_ABBR );    }    }
                              else
                              {    if ( seq_fasta.size( ) == 0 )
                                   {    q.PrintVisual( out, seq[ q.query_id ], qid,
                                             target_fasta[q.target_id], 0,
                                             VISUAL_ABBR );    }
                                   else
                                   {    q.PrintVisual( out, seq_fasta[ q.query_id ],
                                             qid, target_fasta[q.target_id], 0,
                                             VISUAL_ABBR );    }    }    }
                         else
                         {    if ( TARGET_FASTA.size( ) == 0 )
                              {    if ( seq_fasta.size( ) == 0 )
                                   {    q.PrintVisual( out, seq[ q.query_id ], qid,
                                             qual_g[q.target_id],
                                             t, start, VISUAL_ABBR );    }
                                   else
                                   {    q.PrintVisual( out, seq_fasta[ q.query_id ],
                                             qid, qual_g[q.target_id],
                                             t, start, VISUAL_ABBR );    }    }
                              else
                              {    if ( seq_fasta.size( ) == 0 )
                                   {    q.PrintVisual( out, seq[ q.query_id ], qid,
                                             qual_g[q.target_id],
                                             target_fasta[q.target_id], 
                                             0, VISUAL_ABBR );    }
                                   else
                                   {    q.PrintVisual( out, seq_fasta[ q.query_id ],
                                             qid, qual_g[q.target_id],
                                             target_fasta[q.target_id], 
                                             0, VISUAL_ABBR );    }    }    }    }
                    if ( DUMP_WINDOW > 0 )
                    {    seq[id].Print( out, "query_" + ToString(to_query_id(id)) );
                         unsigned int begin = look.ContigStart( q.target_id );
                         static basevector w;
                         unsigned int start = 0;
                         if ( (int) DUMP_WINDOW < q.a.pos2( ) )
                              start = q.a.pos2( ) - DUMP_WINDOW;
                         unsigned int stop = q.a.Pos2( ) + DUMP_WINDOW;
                         if ( !look.CanFetchBasesFromDisk( begin + start,
                              begin + stop ) )
                         {    out << "(Bases from " << begin + start << " to "
                                   << begin + stop
                                   << " are unavailable on disk.  Can't dump "
                                   << "target window.  Sorry.)\n";    }
                         else
                         {    w.Setsize( stop - start );
                              look.FetchBasesFromDisk( begin + start,
                                   begin + stop, w );
                              w.Print( out, "genome_" + target_names[ q.target_id ]
                                   + "." + ToString(start) + "-"
                                   + ToString(stop) );    }    }

                    ++print_count;
                    if ( print_count == (int) KEEP_BEST && print_something )
                    {    out << "...\n";
                         break;    }    }

               // If we found an alignment, the read is placed.  If we found a full
               // length alignment, accept alignments of the read.

               if ( print_count > 0 )
               {    placed[id] = True;
                    if ( these_qualifiers[0].FullLength( ) )
                         placed_full_length[id] = True;
                    if ( print_count == 1 ) uniquely_placed[id] = True;
                    if ( !ALIGN_UNALIGNED_BITS || these_qualifiers[0].FullLength( ) )
                         accepted[id] = True;    }

               i = j - 1;    }

          close(fdq);

          if (LIST_UNPLACED_BY_PASS)
          {    out << "\nUNPLACED QUERY SEQUENCES:";
               int count = 0;
               for ( size_t i = 0; i < seq.size( ); i++ )
               {    int query_id = to_query_id(i);
                    if ( BinMember( to_exclude, query_id ) ) continue;
                    if ( seq[i].size( ) == 0 || seq[i].size( ) < MIN_QUERY_LENGTH )
                         continue;
                    if ( !placed[i] )
                    {    out << ( ( count > 0 && count % 12 == 0 ) ? "\n" : " " );
                         ++count;
                         out << query_id;    }    }
               out << "\n";    }

          int nzseqs = 0;
          for ( size_t x = 0; x < seq.size( ); x++ )
          {    if ( seq[x].size( ) > 0 && seq[x].size( ) >= MIN_QUERY_LENGTH ) 
                    ++nzseqs;    }
          if ( nzseqs > 0 && !QUIET )
          {    out << "\nplacement fraction: " << Sum(placed) << "/" << nzseqs
                    << " = " << PERCENT_RATIO( 3, Sum(placed), nzseqs ) << endl;
               out << "\nfull-length placement fraction: "
                    << Sum(placed_full_length) << "/" << nzseqs << " = "
                    << PERCENT_RATIO( 3, Sum(placed_full_length), nzseqs ) << endl;
               out << "\nunique placement fraction: " << Sum(uniquely_placed)
                    << "/" << nzseqs << " = "
                    << PERCENT_RATIO( 3, Sum(uniquely_placed), nzseqs )
                    << endl;    }    }

     // Report global statistics.

     if (GLOBAL_STATS)
     {    Sort(large_insertions), Sort(large_deletions);
          out << "\nGLOBAL STATISTICS:\n";
          out << total_aligned << " bases on query sequences used in alignments\n";
          out << mismatches << " mismatches\n";
          for ( int i = 1; i < small_insertions.isize( ); i++ )
          {    if ( small_insertions[i] > 0 )
               {    out << small_insertions[i] << " insertions into query of size "
                         << i << "\n";    }    }
          for ( int i = 0; i < large_insertions.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < large_insertions.isize( ); j++ )
                    if ( large_insertions[j] != large_insertions[i] ) break;
               out << j - i << " insertions into query of size "
                    << large_insertions[i] << "\n";
               i = j - 1;    }
          for ( int i = 1; i < small_deletions.isize( ); i++ )
          {    if ( small_deletions[i] > 0 )
               {    out << small_deletions[i] << " deletions from query of size "
                         << i << "\n";    }    }
          for ( int i = 0; i < large_deletions.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < large_deletions.isize( ); j++ )
                    if ( large_deletions[j] != large_deletions[i] ) break;
               out << j - i << " deletions from query of size "
                    << large_deletions[i] << "\n";
               i = j - 1;    }    }

     // Report unplaced sequences.

     if (LIST_UNPLACED)
     {    out << "\nUNPLACED QUERY SEQUENCES:";
          int count = 0;
          for ( size_t i = 0; i < seq.size( ); i++ )
          {    int query_id = to_query_id(i);
               if ( BinMember( to_exclude, query_id ) ) continue;
               if ( seq[i].size( ) == 0 || seq[i].size( ) < MIN_QUERY_LENGTH )
                    continue;
               if ( !placed[i] )
               {    out << ( ( count > 0 && count % 12 == 0 ) ? "\n" : " " );
                    ++count;
                    out << query_id;    }    }
          out << "\n";    }
     if ( UNPLACED_FILE != "" )
     {    Ofstream( un, UNPLACED_FILE );
          for ( size_t i = 0; i < seq.size( ); i++ )
          {    int query_id = to_query_id(i);
               if ( BinMember( to_exclude, query_id ) ) continue;
               if ( seq[i].size( ) == 0 || seq[i].size( ) < MIN_QUERY_LENGTH )
                    continue;
               if ( !placed[i] ) un << query_id << "\n";    }    }
     if ( UNPLACED_SEQUENCE_FILE != "" )
     {    Ofstream( un, UNPLACED_SEQUENCE_FILE );
          for ( size_t j = 0; j < seq.size( ); j++ )
          {    if ( accepted[j] ) continue;
               int query_id = to_query_id(j);
               if ( BinMember( to_exclude, query_id ) ) continue;
               static basevector s;
               s = seq[j];
               if ( seq[j].size( ) < MIN_QUERY_LENGTH ) continue;
               if ( s.size( ) < K ) continue;
               if ( seq_mask[j].size( ) == 0 ) s.Print( un, seq_names[j] );
               else
               {    static vec<ho_interval> covered;
                    covered.clear( );
                    for ( int u = 0; u < seq_mask[j].isize( ); u++ )
                         covered.push_back( ho_interval(
                              seq_mask[j][u].first, seq_mask[j][u].second ) );
                    vec< pair<ho_interval, int> > condensed;
                    CondenseIntervals( s.size( ), covered, condensed );
                    for ( int u = 0; u < condensed.isize( ); u++ )
                    {    if ( condensed[u].second != 0 ) continue;
                         const ho_interval& h = condensed[u].first;
                         if ( h.Length( ) < (int) K ) continue;
                         static basevector b;
                         b.SetToSubOf( s, h.Start( ), h.Length( ) );
                         b.Print( un, seq_names[j] + ".bases:"
                              + ToString( h.Start( ) ) + "-"
                              + ToString( h.Stop( ) ) );    }    }    }    }

     // And we're done!

     if ( p_outStrm != 0 ) delete p_outStrm;
     if ( !QUIET ) out << "\n" << Date( ) << ": DONE!" << endl;
     if (PRINT_MEMORY_USAGE) PrintMemUsage(out);
     if ( OUTFILE != "" ) delete &out;    }
