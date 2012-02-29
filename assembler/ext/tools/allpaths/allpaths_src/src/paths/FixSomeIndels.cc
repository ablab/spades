///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FixSomeIndels.  Attempt to correct some indel error in the assembly.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "Fastavector.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "VecTemplate.h"
#include "efasta/EfastaTools.h"
#include "feudal/QualNibbleVec.h"
#include "math/Functions.h"
#include "paths/AssemblyEdit.h"
#include "paths/FixSomeIndelsUtils.h"
#include "paths/ReadLoc.h"
#include "paths/UnipathFixerTools.h"
#include "util/SearchFastb2Core.h"

void Vote( const vec< vec<int> >& ERRS, const int MIN_INS,
     vec< pair<int,double> >& votes, const Bool VERBOSE, ostringstream *p_rout )
{    
     vec< vec<String> > rows;
     for ( size_t i = 0; i < ERRS.size(); i++ )
     {    vec<int> v = ERRS[i];
          UniqueSort(v);
          if ( v.size( ) > 1 ) 
          {    const vec<int>& errVec = ERRS[i];
               size_t n_err = errVec.size( );
               int best = Max(errVec);
               vec<double> SCORE(n_err);
               double total = 0.0;
               for ( size_t j = 0; j < n_err; j++ )
               {    if ( best - errVec[j] <= 10 )
                    {    SCORE[j] = pow( 10.0, double(errVec[j]-best)/10.0 );
                         total += SCORE[j];    }    }
               for ( size_t j = 0; j < n_err; j++ )
               {    if ( best - errVec[j] <= 10 )
                    {    SCORE[j] /= total;
                         votes.push( MIN_INS + (signed)j, SCORE[j] );    }    }
               vec<String> row;
               if (VERBOSE) 
               {    row.push_back("[" + ToString(i) + "]");
                    bool first = true;
                    String scores, bests;
                    for ( size_t j = 0; j < n_err; j++ ) 
                    {    if ( errVec[j] > -1000 )
                         {    if ( !first ) scores += ", "; 
                              first = false;
                              scores += ToString(MIN_INS + (signed)j) 
                                   + ":" + ToString(errVec[j]);    }    }
                    row.push_back(scores);
                    first = True;
                    for ( size_t j = 0; j < n_err; j++ ) 
                    {    if ( best - errVec[j] <= 10 )
                         {    if ( !first ) bests += ",";
                              first = False;
                              bests += ToString(MIN_INS + (signed)j) 
                                   + ":" + ToString(SCORE[j]);    }    }
                    if ( bests.size( ) > 0 ) row.push_back(bests);
                    rows.push_back(row);    }    }    }
     if (VERBOSE) PrintTabular(*p_rout, rows, 3);    }

void SummarizeAndAcceptVotes( vec< pair<int,double> > votes, 
     vec< pair<double,int> >& votes2, vec<int>& accepted,
     const Bool VERBOSE, ostringstream *p_rout )
{
     // Summarize votes.
      
     Sort(votes);
     for (size_t i = 0; i < votes.size(); i++) 
     {    double X = 0.0;
          size_t j;
          for ( j = i; j < votes.size( ); j++ )
          {    if ( votes[j].first != votes[i].first ) break;
               X += votes[j].second;    }
          votes2.push( X, votes[i].first );
          i = j - 1;    }
     ReverseSort(votes2);
     if (VERBOSE) 
     {    *p_rout << "\nmax votes = " 
               << ( votes2.empty( ) ? 0 : votes2[0].first ) << "\n";
          vec< vec<String> > rows;
          vec<String> row;
          row.push_back( "votes", "delta" );
          rows.push_back(row);
          for (size_t i = 0; i < votes2.size(); i++) 
          {    vec<String> row;
               row.push_back(ToString(votes2[i].first), ToString(votes2[i].second));
               rows.push_back(row);    }
          *p_rout << "\n";
          PrintTabular(*p_rout, rows, 2);    }
      
     // Decide on change and announce.
      
     accepted.clear( );
     double vote_multiplier = 2.5;
     for (size_t i = 0; i < votes2.size(); i++) 
     {    if ( vote_multiplier * votes2[i].first < votes2[0].first ) break;
          accepted.push_back( votes2[i].second );    }    }

void FixIndelsInRepeats( const int tig, const FastaVec & contig, 
     const BaseVecVec & bases, const BaseVecVec & jbases, 
     const VecQualNibbleVec& quals, const VecQualNibbleVec& jquals,
     const vec< triple<int64_t,int,int> > & RALIGNS_this,
     const vec< triple<int64_t,int,int> > & JRALIGNS_this,
     int* p_ambiguities, vec< triple<int,int,String> > * p_edits, 
     ostringstream *p_rout, const String & run_dir, const Bool VERBOSE, 
     const String & TIGS, const int MAX_PERIOD, const double MAX_DEL_FRAC, 
     const double MAX_INS_FRAC, const int MIN_COPIES, const int MIN_LEN,
     const Bool USE_QUALS )
{
     const size_t contig_size = contig.size();
     
     // Find simple sequence repeats.

     for ( size_t cpos = 0; cpos < contig_size; cpos++ )
     {    for ( int period = 1; period <= MAX_PERIOD; period++ )
          {    if ( cpos + (unsigned)MIN_COPIES*period > contig_size ) break;
               Bool amb = False;
               for ( int j = 0; j < period; j++ )
               {    char c = contig[cpos+j];
                    if ( c != 'A' && c != 'C' && c != 'G' && c != 'T' )
                    {    amb = True;
                         break;    }    }
               if (amb) continue;
               int len;
               for ( len = 1; cpos + len < contig_size; len++ )
                    if (contig[cpos + len] != contig[cpos + len % period]) break;
               if ( len < MIN_COPIES*period ) continue;
               if ( len < MIN_LEN ) continue;
               double copies = double(len) / double(period);
               if (VERBOSE) 
               {    *p_rout << "\n";
                    FastaVec seqv;
                    seqv.SetToSubOf( contig, cpos, period );
                    String seq = seqv.ToString( );
                    PRINT6_TO(*p_rout, tig, cpos, period, seq, len, copies);    }
               int rep_start = cpos, rep_stop = cpos + len;
               cpos += len - 1;
      
               FastaVec repeat;
               repeat.SetToSubOf(contig, rep_start, rep_stop - rep_start);
               if (VERBOSE) repeat.Print(*p_rout);
      
               // Select reads.
      
	       vecbasevector reads;
	       VecQualNibbleVec readsq;
	       SelectReads( VERBOSE, USE_QUALS, rep_start, rep_stop,
			    RALIGNS_this, JRALIGNS_this, run_dir, TIGS,
			    bases, jbases, quals, jquals, reads, readsq, p_rout );
	       
               // Try expanding or contracting the repeat.
	       
               if (VERBOSE) 
               {    *p_rout << Date() << ": expanding/contracting repeat of period "
                         << period << endl;    }
               vec< vec<int> > ERRS( reads.size( ) );
               int MIN_INS = int( floor( -MAX_DEL_FRAC * copies ) );
               int MAX_INS = int( floor( MAX_INS_FRAC * copies ) );
               for ( int delta = MIN_INS; delta <= MAX_INS; delta++ ) 
               {    if ( len + delta*period < 0 ) 
                    {    for ( size_t i = 0; i < reads.size( ); i++ )
                              ERRS[i].push_back(-1000);
                         continue;   }
                    FastaVec c1, c2, c3;
                    const int wing = 200; // NEED TO GET RID OF THIS!
                    int cstart = Max( 0, rep_start - wing );
                    int cstop = Min( rep_stop + wing, (int) contig_size );
                    c1.SetToSubOf( contig, cstart, rep_start - cstart );
                    c3.SetToSubOf( contig, rep_stop, cstop - rep_stop );
                    c2.resize( len + delta*period );
                    for ( size_t j = 0; j < c2.size( ); j++ )
                         c2[j] = contig[rep_start + (j % period)];
                    FastaVec c = Cat(c1, c2, c3);
                    const double min_core_identity = 0.9;
                    int nc1 = c1.size( ), nc2 = c2.size( );
                    int agree_size = 2;
                    int NC1 = Max( 0, nc1 - agree_size );
                    int NC2 = Min( nc2 + 2*agree_size, 
                         (int) contig.size( ) - NC1 );
                    FindBestHomeForReads( VERBOSE, ToString( delta ),
					  reads, readsq, c, NC1, NC2, 
					  min_core_identity, ERRS, p_rout );    }
          
               // Compile votes.
               
               if (VERBOSE) *p_rout << "voting (insertion:score)" << endl << endl;
               vec< pair<int,double> > votes;
               Vote( ERRS, MIN_INS, votes, VERBOSE, p_rout );

               // Summarize and accept votes.

               vec<int> accepted;
               vec< pair<double,int> > votes2;
               SummarizeAndAcceptVotes( votes, votes2, accepted, VERBOSE, p_rout );
               const double min_votes = 1.8;
               if ( accepted.solo( ) && accepted[0] == 0 ) break;
               if ( accepted.empty( ) 
                    || ( accepted.nonempty( ) && votes2[0].first < min_votes ) )
               {    accepted.clear( );
                    if ( rep_start == 0 || rep_stop == (int) contig_size )
                         accepted.push_back( 0, -1 );
                    else accepted.push_back( 0, -1, 1 );    }

               // Happy check.  This raised ambiguities a lot, may not be worth it.

               const int min_happy_votes = 4;
               if ( votes2.nonempty( ) && votes2[0].first >= min_votes 
                    && votes2[0].first < min_happy_votes )
               {    for ( int i = 0; i < votes2.isize( ); i++ )
                    {    if ( votes2[i].second == -1 || votes2[i].second == 0
                              || votes2[i].second == 1  
                              || Abs( votes2[i].second - votes2[0].second ) <= 1 )
                         {    accepted.push_back( votes2[i].second );    }    }
                    UniqueSort(accepted);    }

               if (VERBOSE) 
               {    *p_rout << "\nRECOMMEND CHANGE:";
                    for (int i = 0; i < accepted.isize(); i++)
                    *p_rout << " " << accepted[i];
                    *p_rout << "\n";    }
               Sort(accepted);

               *p_ambiguities += (accepted.back() - accepted.front()) * period;

               // Record edit.
      
               EFastaVec e_repeat(repeat.ToString());
               int low = accepted.front();
               if (low > 0) 
               {    for (int l = 0; l < low; l++) 
                    {    for (int i = len - period; i < len; i++)
                         e_repeat.push_back(repeat[i]);    }    }
               else if (low < 0) e_repeat.resize(e_repeat.size() + low*period);
      
               // If the form present in the original contig is in accepted,
               // move it to the front.

               if ( Member( accepted, 0 ) )
               {    vec<int> accepted2;
                    accepted2.push_back(0);
                    for ( int l = 0; l < accepted.isize( ); l++ )
                    {    if ( accepted[l] != 0 ) 
                              accepted2.push_back( accepted[l] );    }
                    accepted = accepted2;    }

               for (int i = 0; i < accepted.isize(); i++)
                    accepted[i] -= low;
               if (accepted.size() > 1) 
               {    
                    // Generate {...} ambiguity notation for accepted.

                    e_repeat.push_back('{');
                    for (int a = 0; a < accepted.isize(); a++) 
                    {    for (int l = 0; l < accepted[a]; l++) 
                         {    for (int i = len - period; i < len; i++)
                                   e_repeat.push_back(repeat[i]);    }
                         if (a < accepted.isize() - 1) e_repeat.push_back(',');    }
                    e_repeat.push_back('}');    }
      
               p_edits->push(rep_start, rep_stop, e_repeat);
      
               break;    }    }    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_Int(K);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
	   "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault(HEAD, "extended40.shaved");
     CommandArgument_String_OrDefault(SCAFFOLDS_IN, "linear_scaffolds0.patched");
     CommandArgument_String_OrDefault(SCAFFOLDS_OUT, 
          "linear_scaffolds0.patched.fixed");
     CommandArgument_String_OrDefault(UNIPATH_PATCH_DIR, "unipath_patch");
     CommandArgument_String_OrDefault(POST_PATCH_DIR, "post_patch");

     CommandArgument_String_OrDefault_Doc(TIGS, "", 
          "if unspecified, process all contigs;"
          " otherwise it is one of the following: \n"
          "(a) a list of contig ids (in ParseIntSet format) or \n"
          "(b) the letter s followed by a list of scaffolds or \n"
          "(c) s<scaffold id>.<list of indices of contigs in the scaffold");
     CommandArgument_String_OrDefault_Doc(TARGETS, "", "for assessment, list of "
          "reference targets to use, in ParseIntSet format, to speed up alignment");
     CommandArgument_String_OrDefault(FRAG_READS, "frag_reads_filt");
     CommandArgument_String_OrDefault(JUMP_READS, "jump_reads_filt_cpd");
     CommandArgument_Int_OrDefault(MAX_PERIOD, 8);
     CommandArgument_Double_OrDefault(MAX_DEL_FRAC, 0.55);
     CommandArgument_Double_OrDefault(MAX_INS_FRAC, 0.56);
     CommandArgument_Int_OrDefault(MIN_COPIES, 7);
     CommandArgument_Bool_OrDefault(USE_QUALS, True);
     CommandArgument_Int_OrDefault(MIN_LEN, 7);
     CommandArgument_Bool_OrDefault(SHOW_FLAKY, False);
     CommandArgument_Bool_OrDefault(SHOW_FIX_INDELS, False);
     CommandArgument_Bool_OrDefault(WRITE, True);
     EndCommandArguments;

     // Begin.

     double clock = WallClockTime();

     // Thread control

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Define directories, etc.

     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     cout << Date() << ": " << run_dir << endl;
     Bool CHECKPOINT = True;
     String ch_head = sub_dir + "/" + POST_PATCH_DIR + "/PostPatcher." + SCAFFOLDS_IN + ".";
     Mkpath( sub_dir + "/" + POST_PATCH_DIR );

     // Set up read locations.  Load the reads.

     read_locs_on_disk locs_file( sub_dir + "/" + SCAFFOLDS_IN, run_dir );
     vecbasevector frag_reads( run_dir + "/" + FRAG_READS + ".fastb" );
     VecQualNibbleVec frag_quals;
     LoadQualNibbleVec( run_dir + "/" + FRAG_READS + ".qualb", &frag_quals );

     // Note some files that are needed.

     String unibases_file = run_dir + "/" + HEAD + ".unibases.k" + ToString(K);
     String TIGS_file = ch_head + "TIGS";
     String tigsa_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fasta";
     String jreads_file = run_dir + "/" + JUMP_READS + ".fastb";

     String out_head = sub_dir + "/" + SCAFFOLDS_OUT;

     String scaffolds_tig_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.vecfasta";

     String psort_file = ch_head + "JRALIGNS_PSORT";
     ForceAssert(IsRegularFile(psort_file));

     // Get total number of contigs.
   
     const size_t n_tigs = MastervecFileObjectCount(scaffolds_tig_file);

     // Parse arguments.

     vec<int> tigs_in;
     {    if (TIGS == "") 
          {    for (size_t j = 0; j < n_tigs; j++)
                    tigs_in.push_back(j);    }
          else if (TIGS.Contains("s", 0)) 
          {    TIGS = TIGS.After("s");
               vec<superb> scaffolds;
               ReadSuperbs(sub_dir + "/" + SCAFFOLDS_IN + ".superb", scaffolds);
               if (TIGS.Contains(".")) 
               {    int scaffold = TIGS.Before(".").Int();
                    ForceAssertLt(scaffold, scaffolds.isize());
                    vec<int> spos;
                    ParseIntSet(TIGS.After("."), spos);
                    for (int j = 0; j < spos.isize(); j++)
                         tigs_in.push_back(scaffolds[scaffold].Tig(spos[j]));    }
               else 
               {    vec<int> s;
                    ParseIntSet(TIGS, s);
                    for (int i = 0; i < s.isize(); i++) 
                    {    int scaffold = s[i];
                         ForceAssertLt(scaffold, scaffolds.isize());
                         for (int j = 0; j < scaffolds[scaffold].Ntigs(); j++)
                              tigs_in.push_back(scaffolds[scaffold].Tig(j));
                              }    }    }
          else ParseIntSet(TIGS, tigs_in);    }

     // Create index "ids" that corresponds to putting the contigs in descending
     // order by size.
 
     const vec<int>& tigs = tigs_in;
     const size_t n_used_tigs = tigs.size( );
     vec<int> len(n_used_tigs), ids( n_used_tigs, vec<int>::IDENTITY );
     for ( size_t i = 0; i < n_used_tigs; i++ )
     {    FastaVecVec contigs;
          contigs.ReadOne(scaffolds_tig_file, tigs[i]);
          len[i] = contigs[0].size( );    }
     ReverseSortSync( len, ids );

     vec<int> left_trim( n_tigs, 0 ), right_trim( n_tigs, 0 );

     // Set up to track original and expanded contigs (for evaluation purposes).
  
     vec< vec<FastaVec> > contigs_orig(n_used_tigs);
     vec<EFastaVec> econtigs_new(n_used_tigs);
     // Collect all the edits and save in external file at the end, so that we can 
     // combine edits from different error-correction pograms (currently FixSomeIndels 
     // and FixAssemblyBaseErrors) without re-alignment
     vec< vec< triple< int, int, String > > > all_edits(n_used_tigs); // (start, stop. replacement)
  
      // Create position-sorted versions of read alignment files and indices to them.

     {
  
      vec<size_t> JRPS_START;
      BinaryRead3(ch_head + "JRPS_START", JRPS_START);
  
      // Load the reads.
  
      BaseVecVec jbases;
      VecQualNibbleVec jquals;
      if (TIGS == "") {
        cout << Date() << ": loading jump reads" << endl;
        jbases.ReadAll(run_dir + "/" + JUMP_READS + ".fastb");
        cout << Date() << ": reads loaded, mem = "
             << ToStringAddCommas(MemUsageBytes()) << endl;
        if (USE_QUALS)
	  LoadQualNibbleVec( run_dir + "/" + JUMP_READS + ".qualb", &jquals );
      }

      int processed = 0;
      cout << Date() << ": " << n_used_tigs << " contigs to process" << endl;
      int tigs_per_dot = (TIGS == "" ? 1000 : 10);
      cout << Date() << ": one dot per " << tigs_per_dot 
           << " contigs processed:" << endl;

     size_t total_bases = 0;
     vec<int> ambiguities(n_used_tigs, 0);
     vec<String> report(n_used_tigs);
     Bool SHOW = ( SHOW_FIX_INDELS || SHOW_FLAKY );
  
     #pragma omp parallel for schedule(dynamic, 1)
     for (size_t it = 0; it < n_used_tigs; it++) 
     {    int itx = ids[it];
          int tig = tigs[itx];
          ostringstream rout;
          if (SHOW) rout << "\n" << Date() << ": tig = " << tig << endl;
    
          // Load the contig.

          FastaVecVec contigs;
          vec< triple<int64_t, int, int> > JRALIGNS_this;
          #pragma omp critical
          {
            contigs.ReadOne(scaffolds_tig_file, tig);
            contigs_orig[itx].push_back(contigs[0]);

            // Get the read alignments for this contig.

            if (SHOW) 
              rout << Date( ) << ": fetching jump reads for this contig" << endl;
            BinaryReadRange3(psort_file,JRPS_START[tig], JRPS_START[tig+1], JRALIGNS_this);
          }
          FastaVec& contig = contigs[0];

          // Trim unreliable bases off ends of contigs.  At present this uses
          // only the evidence provided by fragment reads.  If a base in the first
          // or last 60 bases of a contig is less than 61% supported by the reads,
          // we trim off that base and all the bases from it to the end.
          // Also display reads around flaky regions.

          vec<read_loc> locs;
          #pragma omp critical
          {    locs_file.LoadContig( tig, locs );    }
          vec<Bool> to_remove( locs.size( ), False );
          for ( size_t i = 0; i < locs.size( ); i++ )
               if ( locs[i].ReadClass( ) != 0 ) to_remove[i] = True;
          EraseIf( locs, to_remove );
          basevector TIG = contig.ToBasevector( );
          int n = TIG.size( );
          vecbasevector jump_reads, long_jump_reads;
          vec<dumbcall> calls;
          Pileup( TIG, locs, frag_reads, jump_reads, long_jump_reads, calls );
          const double min_ref_frac = 0.61;
          if (SHOW_FLAKY)
	    ShowFlakyRegions( n, tig, min_ref_frac, TIG, frag_reads, calls, locs, rout );
          const int head = 60;
          int left, right; // number of bases to trim from left and right
          for ( left = Min( head, n - 1 ); left >= 0; left-- )
          {    double ref_frac = double( calls[left].Count( TIG[left] ) ) 
                    / double( calls[left].CountAll( ) );
               if ( ref_frac < min_ref_frac ) break;    }
          for ( right = Min( head, n - 1 ); right >= 0; right-- )
          {    double ref_frac = double( calls[n-right-1].Count( TIG[n-right-1] ) ) 
                    / double( calls[n-right-1].CountAll( ) );
               if ( ref_frac < min_ref_frac ) break;    }
          left++, right++;
          PRINT3_TO( rout, tig, left, right );

          // It's possible that a contig is trimmed to nothing.  We don't actually
          // do this, becaused it could result in zero-length scaffolds, which
          // would break subsequent code.

          if ( left + right >= (int) contig.size( ) ) left = right = 0;

          // Now record trims and edit the contig.

          left_trim[tig] = left, right_trim[tig] = right;
          contig.SetToSubOf( contig, left, (int) contig.size( ) - left - right );
          
	  // Set up to track edits.

	  //vec< triple<int, int, String> > edits; // (start, stop, replacement)
	  vec< triple<int, int, String> >& edits = all_edits[itx]; // (start, stop, replacement)

	  // Fix indels in simple sequence repeats.

          if (SHOW_FIX_INDELS) rout << Date( ) << ": FixIndelsInRepeats:" << endl;
          vec< triple<int64_t,int,int> > frag_RALIGNS_this;
          for ( int j = 0; j < locs.isize( ); j++ )
          {    const read_loc& rl = locs[j];
               frag_RALIGNS_this.push( rl.ReadId( ), rl.ContigId( ),
                    ( rl.Fw( ) ? rl.Start( ) - left 
                    : -rl.Start( ) - 1 + left ) );     }
          for ( int j = 0; j < JRALIGNS_this.isize( ); j++ )
          {    int& p = JRALIGNS_this[j].third;
               if ( p >= 0 ) p -= left;
               else p += left;    }
          // ARGH: NOTE THAT THIS DOESN'T CORRECTLY HANDLE SIGNS....
     
          FixIndelsInRepeats( tig, contig, frag_reads, jbases, frag_quals, jquals,
               frag_RALIGNS_this, JRALIGNS_this, &ambiguities[itx], &edits, &rout,
               run_dir, SHOW_FIX_INDELS, TIGS, MAX_PERIOD, MAX_DEL_FRAC, 
	       MAX_INS_FRAC, MIN_COPIES, MIN_LEN, USE_QUALS );

          size_t n_edits = edits.size();

          // Sanity check.  Note that we can't leave this here.  Instead need
          // to change the previous code so that it's guaranteed.

          for (size_t ie = 0, je = 1; je < n_edits; ie++, je++)
               ForceAssert(edits[ie].second <= edits[je].first);

          // Make edits.
          size_t contig_size = contig.size(), ie = 0;
	  if ( WRITE ) // Will not make the change until WRITE switch is on
	  {
	    for ( size_t i = 0; i < contig_size; i++ ) 
	    { 
	      if ( ie < n_edits && i == (unsigned)edits[ie].first )
	      {    
		for (size_t j = 0; j < edits[ie].third.size(); j++)
		  econtigs_new[itx].push_back( edits[ie].third[j] );
		i += edits[ie].second - edits[ie].first - 1;
		ie++;    
	      }
	      else  
		if (  contig[i] == 'A' || contig[i] == 'C' || contig[i] == 'G' || contig[i] == 'T' )
		  econtigs_new[itx].push_back( contig[i] );    
		else
		  econtigs_new[itx].append( ExpandAmbCode(contig[i]) );  
	    }
	  }
	  else // Make no edits, all the edit positions are adjust for the original un-trimmed contig
	  {
	    for( size_t k = 0; k < edits.size(); k++ )
	    {
	      edits[k].first += left;
	      edits[k].second += left;
	    }
	    if ( left > 0 )
	      edits.push( 0, left, "" );
	    if ( right > 0 )
	      edits.push( contig.size() + left, contig.size() + left +  right,  "");
	  }
          if (SHOW) report[itx] = rout.str();

          // Record progress.

          #pragma omp critical
          {    total_bases += contig_size;
               DotMod(cout, processed, tigs_per_dot);
               processed++;    }    }
      cout << endl;

     // Report results.
    
     cout << Date() << ": done processing contigs" << endl;
     if (SHOW) 
     {    for (size_t it = 0; it < n_used_tigs; it++)
               cout << report[it];    }
     cout << "bases trimmed = " << Sum(left_trim) + Sum(right_trim) << endl;
     cout << "added ambiguities = " << ToStringAddCommas(Sum(ambiguities)) << endl;
     cout << "total bases = " << ToStringAddCommas(total_bases) << endl;
     
    } // unneeded data gets destroyed here

     // Write output

     if ( ! WRITE ) {
       String outputFile =  out_head + ".edits";
       cout << "Assembly will not be saved (running with WRITE=False)\n" << endl;
       cout << "Output suggested edits in file " << outputFile << endl;
       vec<assembly_edit> edits;
       for (size_t it = 0; it < n_used_tigs; it++)
       {    for( size_t k = 0; k < all_edits[it].size(); k++ )
	    {    vec<basevector> reps;
                 efasta(all_edits[it][k].third).ExpandTo(reps);
                 edits.push( assembly_edit::INTERNAL, tigs[it], 
                      all_edits[it][k].first, 
                      all_edits[it][k].second, reps );    }    }
       BinaryWriter::writeFile( outputFile.c_str( ), edits );
       return 0;
     }

     cout << "\n" << Date( ) << ": Writing fixed efasta files to : " 
          << out_head << ".*" << endl;
     vec<FastaVec> flattened_fasta(econtigs_new.size());
     vec<FastaVec> flattened_fasta_max(econtigs_new.size());
     vecbasevector flattened_fastb(econtigs_new.size());
     {    Ofstream(out_e, out_head + ".contigs.efasta");
          Ofstream(out_a, out_head + ".contigs.fasta");
	  Ofstream(out_m, out_head + ".contigs.max.fasta");
          for (size_t it = 0; it < n_used_tigs; it++) 
          {    econtigs_new[it].FlattenTo(flattened_fasta[it]);
	       econtigs_new[it].FlattenNMaxTo(flattened_fasta_max[it]);
               econtigs_new[it].FlattenTo(flattened_fastb[it]);
               econtigs_new[it].Print(out_e, ToString(it));
               flattened_fasta[it].Print(out_a, ToString(it));   
	       flattened_fasta_max[it].Print(out_m, ToString(it));  }
          flattened_fastb.WriteAll( out_head + ".contigs.fastb" );     }
     vecfastavector flattened_fasta2(econtigs_new.size());
     for (size_t it = 0; it < n_used_tigs; it++)
          flattened_fasta2.push_back_reserve( flattened_fasta[it] );
     flattened_fasta2.WriteAll( out_head + ".contigs.vecfasta" );
     vec<superb> scaffolds;
     ReadSuperbs(sub_dir + "/" + SCAFFOLDS_IN + ".superb", scaffolds);
     for (int ii=0; ii<scaffolds.isize( ); ii++) 
     {    superb& S = scaffolds[ii];
          for ( int jj=0; jj < S.Ntigs( ); jj++ ) 
          {    int tig_id = S.Tig( jj );
               int tig_len = flattened_fastb[tig_id].size( );
               S.SetLen( jj, tig_len );    
               if ( jj > 0 ) S.SetGap( jj-1, S.Gap(jj-1) + left_trim[tig_id] );
               if ( jj < S.Ngaps( ) - 1 ) 
                    S.SetGap( jj, S.Gap(jj) + right_trim[tig_id] );    }    }
     WriteScaffoldedEFasta(out_head + ".assembly.efasta", econtigs_new, scaffolds);
     WriteScaffoldedFasta(out_head + ".assembly.fasta", flattened_fasta, scaffolds);
     WriteScaffoldedFasta(out_head + ".assembly.max.fasta", flattened_fasta_max, scaffolds);
     WriteSuperbs( out_head + ".superb", scaffolds );
     WriteSummary( out_head + ".summary", scaffolds );
     
     // Done.

     cout << "time used = " << TimeSince(clock) << endl;    }
