///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FixAssemblyBaseErrors.  Attempt to correct some indel error in the
// assembly (code copied from FixSomeIndels).

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "Fastavector.h"
#include "MainTools.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "VecTemplate.h"
#include "efasta/EfastaTools.h"
#include "feudal/QualNibbleVec.h"
#include "math/Functions.h"
#include "paths/AssemblyEdit.h"
#include "paths/CAltFasta.h"
#include "paths/FixSomeIndelsUtils.h"
#include "paths/ReadLoc.h"


// Given the base call frequency, calculate the fraction of the 
// calls same as the reference base
double AgreeFraction(const dumbcall& call, int ref_base)
{
  int ncalls = 0, ncalls_agree = 0;
  for (int k = 0; k < 6; k++)
  {   
    ncalls += call.base[k];
    if ( k == ref_base ) ncalls_agree += call.base[k];    
  }
  return double(ncalls_agree)/double(ncalls);
}

// Given the base call frequencies, return the possible base
// choices;
vec<fastavector> AltBases(const dumbcall& call, const int ref_base)
{
  ForceAssert( ref_base >= 0 && ref_base < 4 );
  double cutoff = 0.20; // sy_fix_later
  vec< fastavector > results;
  results.push_back( fastavector( String(as_base(ref_base)) ) );
  int ncalls = call.CountAll();
  for (int k = 0; k < 6; k++)
  {   
    if ( double( call.Count(k) ) / ncalls  < cutoff ) continue;
    if ( k == ref_base ) continue;
    if ( k == 4 ) // 'D', get empty fasta string
    {
      results.push_back( fastavector("") );
    }
    else if ( k == 5 ) // insertion, add all 4 posibilities choices {A,C,G,T}<original>
    {
      // sy_fix_later
      //for ( int kk = 0; kk < 4; k++ )
      //{
      //  fastavector bases;
      //  bases.push_back( as_base(kk) );
      //  bases.push_back( as_base(ref_base) );
      //  results.push( bases );
      //}
    }
    else
    {
      results.push_back( fastavector( String(as_base(k)) ) );
    }
  }
  return results;
}


void GetFlakySites(const fastavector& contig, const vec<read_loc> locs, const int tig,
    const vecbasevector& frag_reads, const vecbasevector& jump_reads, const vecbasevector& long_jump_reads,
   vec<int>& flakySites, vec<dumbcall>& calls, ostream& rout = std::cout, bool VERBOSE = false )
{
  flakySites.clear();
  calls.clear();
  basevector TIG = contig.ToBasevector( );
  Pileup( TIG, locs, frag_reads, jump_reads, long_jump_reads, calls );
  const double min_ref_frac = 0.60;  // sy_fix_later
  // determine the percentage of the calls agree with the reference
  for (int j = 0; j < TIG.isize( ); j++)
  {
    // no ambigous sites
    if ( contig[j] != 'A' &&  contig[j] != 'C' &&  contig[j] != 'G' &&  contig[j] != 'T' )
      continue;
    if (AgreeFraction(calls[j], TIG[j]) < min_ref_frac) 
    {
      if ( AltBases( calls[j], TIG[j] ).size()  > 1 )
      {
	flakySites.push(j);
	if ( VERBOSE ) 
	  rout << "add site " << j << " " << contig[j] << " " << AgreeFraction(calls[j], TIG[j]) << " " << calls[j].CountAll() << endl;
      }
      else
	if ( VERBOSE ) 
	  rout << "    skip " << j << " " << contig[j] << " " << AgreeFraction(calls[j], TIG[j]) << " " << calls[j].CountAll() << endl;
    }
  }
  if ( VERBOSE ) rout << "Total flaky sites " << flakySites.size() << endl;
}

int GetAltFasta(const fastavector& contig, const int tig, const vec<int>& flakySites, const vec<dumbcall>& calls,
   vec<CAltFasta>& results, ostream& rout = std::cout, bool VERBOSE = false )
{
  basevector TIG = contig.ToBasevector( );
  int tlen = TIG.size();
  // find the segments that need to be corrected, which satisfy the following rules
  // 1. There are more than 20 flaky sites on both sides ( radius = 15)
  // 2. The total length of the segments is less than 100
  // 3. The number of ambiguous sites are less than 3
  int radius = 20, MaxFlakySites = 8, MaxSegLength = 100; // sy_fix_later
  int left = -1, right = -1;
  int left_wing = -1, right_wing = -1; // how much safe zone on both sides
  int n_esites = 0; // total number of sites to be evaluated
  for (int i = 0; i < flakySites.isize(); i++)
  {
    if ( VERBOSE ) rout << "site " << i <<" "<<flakySites[i] << endl;
    int sep_left_i = ( i > 0 ) ?  flakySites[i] - flakySites[i-1] : flakySites[i] - 0;
    if ( sep_left_i > radius ) { left = i; left_wing = sep_left_i; }
    right_wing = ( i <= flakySites.isize() ) ? flakySites[i+1] - flakySites[i] : tlen -  flakySites[i];
    if ( right_wing > radius ) 
    {
      right = i;
      if ( left < 0 ) continue;
      int start = flakySites[left] - 0;
      int end =  flakySites[right] + 0 + 1;  // [start, end) is the region needs to be corrected
      if ( right - left + 1 > MaxFlakySites || end - start >= MaxSegLength)
       	continue;
      // prepare the output
      fastavector seg = contig;
      vec<vec<fastavector> > candBases;
      vec<int> pickSites;
      for(int k = left; k <= right; k++)
      {
	int site = flakySites[k];
	candBases.push_back( AltBases( calls[site], TIG[site] ) );
	pickSites.push_back( site - start );
      }
      n_esites += pickSites.size();
      if ( VERBOSE ) 
      {
	rout << "    " << start << " " << end << " " << pickSites.size() << endl;
	for(size_t k = 0; k< pickSites.size(); k++ )
	{ 
	  int site = pickSites[k] + start; 
	  rout <<  "     " << site << ":" << contig[site] << "(" << as_base(TIG[site]) << ") -> " ;
	  for ( size_t kk = 0; kk < candBases[k].size(); kk++ )
	  {
	    String str =  candBases[k][kk].ToString();
	    if ( str == "" ) str = "D";
	    rout << str << " ";
	  }
	  rout << endl;
	}
      }
      int select_left_wing = Min( left_wing, 200);
      int select_right_wing = Min( right_wing, 200);
      results.push( CAltFasta( tig, start, end, contig, pickSites, candBases, select_left_wing , select_right_wing) );
      left = -1;
    }
  }
  return n_esites;
} 


int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_Int(K);
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
     CommandArgument_String_OrDefault(JUMP_READS, "jump_reads_filt_cpd");
     CommandArgument_Bool_OrDefault(USE_QUALS, True);
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
	   "Number of threads to use (use all available processors if set to 0)");
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
     vecbasevector frag_reads( run_dir + "/frag_reads_filt.fastb" );
     VecQualNibbleVec frag_quals;
     if (USE_QUALS)
       LoadQualNibbleVec( run_dir + "/frag_reads_filt.qualb", &frag_quals );

     // Note some files that are needed.

     String unibases_file = run_dir + "/" + HEAD + ".unibases.k" + ToString(K);
     String TIGS_file = ch_head + "TIGS";
     String tigsa_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fasta";
     String jreads_file = run_dir + "/" + JUMP_READS + ".fastb";

     String out_head = sub_dir + "/" + SCAFFOLDS_OUT;

     String scaffolds_tig_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.vecfasta";

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
      {
        cout << Date() << ": loading jump reads" << endl;
        jbases.ReadAll(run_dir + "/jump_reads_filt_cpd.fastb");
        cout << Date() << ": reads loaded, mem = "
             << ToStringAddCommas(MemUsageBytes()) << endl;
        if (USE_QUALS)
	  LoadQualNibbleVec( run_dir + "/jump_reads_filt_cpd.qualb", &jquals );
      }

      int processed = 0;
      cout << Date() << ": " << n_used_tigs << " contigs to process" << endl;
      int tigs_per_dot = (TIGS == "" ? 1000 : 10);
      cout << Date() << ": one dot per " << tigs_per_dot 
           << " contigs processed:" << endl;

     size_t total_bases = 0;
     int total_fsites = 0, total_esite = 0, total_csites = 0;

     vec<int> ambiguities(n_used_tigs, 0);
     vec<String> report(n_used_tigs);
     Bool SHOW =  VERBOSE;
  
     #pragma omp parallel for schedule(dynamic, 1)   
     for (size_t it = 0; it < n_used_tigs; it++) 
     {    
       int itx = ids[it];
       int tig = tigs[itx];
       //ostream& rout = cout;
       ostringstream rout;
       if (SHOW) rout << "\n" << Date() << ": tig = " << tig << endl;
       
       int n_fsites = 0, n_csites = 0, n_esites =0;
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
	 BinaryReadRange3(ch_head + "JRALIGNS_PSORT", 
	     JRPS_START[tig], JRPS_START[tig+1], JRALIGNS_this);
       }
       if (SHOW) rout << Date( ) << ": End fetching jump reads " << endl;
       FastaVec& contig = contigs[0];

       vec<read_loc> locs;
       #pragma omp critical
       {    locs_file.LoadContig( tig, locs );    }
       vec<Bool> to_remove( locs.size( ), False );
       for ( size_t i = 0; i < locs.size( ); i++ )
         //if ( locs[i].ReadClass( ) != 0 ) to_remove[i] = True;
         if ( locs[i].ReadClass( ) == 2 ) to_remove[i] = True;
       EraseIf( locs, to_remove );
       vecbasevector jump_reads, long_jump_reads;
       if (SHOW) rout << Date( ) << ": End fetching frag reads " << endl;

       // Find all flaky sites and the pile up of reads
       vec<CAltFasta> altFastas;
       {
	 vec<int> flakySites;
	 vec<dumbcall> calls;
	 //GetFlakySites(contig, locs, tig, frag_reads, jump_reads, long_jump_reads, flakySites, calls, rout, VERBOSE);
	 GetFlakySites(contig, locs, tig, frag_reads, jbases, long_jump_reads, flakySites, calls, rout, VERBOSE);
	 // not prepare the segments with alternative bases for re-evaluation
	 n_fsites = flakySites.size();
	 n_esites = GetAltFasta(contig, tig, flakySites, calls, altFastas, rout, VERBOSE);
       }
       // Collect reads in the region
       if ( SHOW ) rout << "alternative fasta segments = " << altFastas.size() << endl;
       if ( SHOW ) rout << "alternative fasta sites = " << n_esites << endl;
       // aprepare the lignments of fragment reads
       vec< triple<int64_t,int,int> > frag_RALIGNS_this;
       for ( int j = 0; j < locs.isize( ); j++ )
       {    
	 const read_loc& rl = locs[j];
	 if ( rl.ReadClass( ) == 0 )
	   frag_RALIGNS_this.push( rl.ReadId( ), rl.ContigId( ),
	       ( rl.Fw( ) ? rl.Start( ) : -rl.Start( ) - 1 ) );     
	 // sy_fix_later
	 if (  rl.ReadClass( ) == 1 )
	   JRALIGNS_this.push( rl.ReadId( ), rl.ContigId( ),
	       ( rl.Fw( ) ? rl.Start( ) : -rl.Start( ) - 1 ) );     
       }
       //Sort( JRALIGNS_this );
       UniqueSort(JRALIGNS_this);
       Sort(JRALIGNS_this, cmp_pos);
       // All edits 
       //vec< triple<int, int, String> > edits; // (start, stop, replacement)
       vec< triple<int, int, String> >& edits = all_edits[itx]; // (start, stop, replacement)
       for( size_t ia = 0; ia < altFastas.size(); ia++ )
       {
	 altFastas[ia].SetOstringstream( rout );
	 if ( SHOW ) 
	   altFastas[ia].Print( rout );
	 //altFastas[ia].EvalAlternatives( false, USE_QUALS, TIGS, run_dir, frag_reads, jbases, frag_quals, jquals, frag_RALIGNS_this, JRALIGNS_this );
	 altFastas[ia].EvalAlternatives( VERBOSE, frag_reads, jbases, frag_quals, jquals, frag_RALIGNS_this, JRALIGNS_this );
	 altFastas[ia].Vote( VERBOSE );
	 altFastas[ia].SummarizeAndAcceptVotes( VERBOSE );
	 altFastas[ia].AcceptedEdits( edits, VERBOSE );
	 //if ( SHOW ) 
	 //  altFastas[ia].Print( rout );
       }
       n_csites = edits.size();

       Sort(edits);
       size_t n_edits = edits.size();
       if ( VERBOSE ) 
       {
	 rout << Date() << " Collecting edits for tig" << tig << endl;
	 rout << "n_edits= " << n_edits << endl;
	 for ( size_t i = 0; i < edits.size(); i++ )
	 {
	   rout << ">tig" << tig << "-" << i <<" " << edits[i].first << "," << edits[i].second << " " 
	     << contig[edits[i].first] << " -> " << edits[i].third << endl;
	   fastavector temp;
	   temp.SetToSubOf( contig, Max( 0, edits[i].first - 100 ), Min( 200u, contig.size() - edits[i].first + 100 ) );
	   rout << temp.ToString() << endl;
	   temp.ReverseComplement();
	   rout << temp.ToString() << endl;
	 }
       }
       // Will not make the change until WRITE switch is on
       size_t contig_size = contig.size(), ie = 0;
       if ( WRITE ) 
       { // make the edits
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
       // collect report and statistics
       if (SHOW) report[itx] = rout.str();
#pragma omp critical
       {    
	 total_fsites += n_fsites;
	 total_esite += n_esites;
	 total_csites += n_csites;
	 total_bases += contig_size;
	 DotMod(cout, processed, tigs_per_dot);
	 processed++;    
       }    
     } // end for each contig
     cout << endl;
     // Report results.
     cout << Date() << ": done processing contigs" << endl;
     if (SHOW) 
     {    for (size_t it = 0; it < n_used_tigs; it++)
               cout << report[it];    }
     cout << "added ambiguities = " << ToStringAddCommas(Sum(ambiguities)) << endl;
     cout << "total bases = " << ToStringAddCommas(total_bases) << endl;
     cout << "total flaky sites = " << ToStringAddCommas(total_fsites) << endl;
     cout << "total evaluated sites = " << ToStringAddCommas(total_esite) << endl;
     cout << "total changed sites = " << ToStringAddCommas(total_csites) << endl;
     
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
          {    
	       econtigs_new[it].FlattenTo(flattened_fasta[it]);
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
