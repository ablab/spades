///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
   PathsToLocs
   
   Given some kmer paths, find all their possible placements on the genome.  
   None of the kmer paths are allowed to have gaps.

   We assume that a consistent numbering for the kmer paths and the genome has
   been provided.  If this is not the case, and genome.fastb exists, use the 
   RENUMBER option.  Otherwise, you'll get garbage.
*/


#include "Basevector.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "graph/Digraph.h"
#include "graphics/Whiteboard.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"
#include "paths/simulation/GenomePlacement.h"
#include "paths/PdfEntry.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN); 
     CommandArgument_String_OrDefault(PATHS, ""); 
     CommandArgument_Int(K);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
       "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault(SAVE_PLACEMENTS_TO, "");
     CommandArgument_String_OrDefault(LOAD_PLACEMENTS_FROM, "");
     CommandArgument_Bool_OrDefault(UNIPATH, False);
     CommandArgument_Bool_OrDefault(FW_ONLY, False);
     CommandArgument_String_OrDefault(COPY_PDF_FILE, "");
     CommandArgument_Int_OrDefault(MATCH_LEN_PERCENT, (UNIPATH ? 75 : 100) );
     CommandArgument_Bool_OrDefault(SHOW_PLACEMENTS, !UNIPATH);
     CommandArgument_Bool_OrDefault(SHOW_FW_ONLY, True);
     CommandArgument_String_OrDefault(SHOW_THESE_CONTIGS_ONLY, "");
     CommandArgument_Int_OrDefault(SHOW_WINDOW_BEGIN, -1);
     CommandArgument_Int_OrDefault(SHOW_WINDOW_END, -1);
     CommandArgument_Bool_OrDefault(KMER_COORDS, False);
     CommandArgument_Int_OrDefault(MIN_KMERS, 1);
     CommandArgument_Int_OrDefault(MAX_KMERS, 1000000000);
     CommandArgument_Bool_OrDefault(UNIPATH_PLOT, False);
     CommandArgument_Bool_OrDefault(PRINT_ADJACENCY_GRAPH, True);
     CommandArgument_Int_OrDefault_Doc(PRINT_UNPLACED, 1000000000,
          "print unplaced unipaths having at least the given size");
     CommandArgument_Int_OrDefault(CUTOFF, 0);
     CommandArgument_Bool_OrDefault(RENUMBER, False);
     CommandArgument_Bool_OrDefault(RC_TOO, False);
     CommandArgument_String_OrDefault(READS, "reads");
     CommandArgument_String_OrDefault(UNIPATHS, "unipaths");
     CommandArgument_String_OrDefault(UNIPATHS_FASTB, "");
     CommandArgument_Bool_OrDefault(PATHS_ABSOLUTE, False); 
     CommandArgument_String_OrDefault_Doc(ALT_DATA_DIR, "",
          "alternate data_dir, allowing one to substitute a different reference");
     CommandArgument_String_OrDefault(GENOME, "genome");
     CommandArgument_Bool_OrDefault(PRINT_GAPS, True); 
     CommandArgument_Bool_OrDefault(STATS, True); 
     CommandArgument_String_OrDefault(CN_SUFFIX, "");
     CommandArgument_Bool_OrDefault(COPY_NUMBER_UNKNOWN, False);
     CommandArgument_String_OrDefault_Doc(CN_TIGS, "",
          "reference contigs to use for copy number evaluation");
     CommandArgument_Int_OrDefault(TIG_OFFSET, 0);
     CommandArgument_Int_OrDefault(POS_OFFSET, 0);
     CommandArgument_Bool_OrDefault(ARCHIVE, False);
     EndCommandArguments;

     // Thread control
     
     NUM_THREADS = configNumThreads(NUM_THREADS);
 
     // Define directories.

     String datadir = PRE + "/" + DATA;
     if ( ALT_DATA_DIR != "" ) datadir = ALT_DATA_DIR;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String tmp_dir = run_dir + "/tmp";
     Mkdir777(tmp_dir);

     // Define copy number file.

     String copy_pdf_file( COPY_PDF_FILE );
     if( copy_pdf_file.empty() )
     {    if ( CN_SUFFIX != "" ) CN_SUFFIX = "." + CN_SUFFIX;
          copy_pdf_file = run_dir + "/" + READS + "." + UNIPATHS 
               + ".predicted_count.k" + ToString(K) + CN_SUFFIX;    }

     // Sanity check arguments.

     if( MATCH_LEN_PERCENT < 50 ) FatalErr("MATCH_LEN_PERCENT must be at least 50%");
     if( SAVE_PLACEMENTS_TO.nonempty() && LOAD_PLACEMENTS_FROM.nonempty() )
     {    FatalErr( "Specify at most one of SAVE_PLACEMENTS_TO or "
               "LOAD_PLACEMENTS_FROM filenames" );    }
     if( PATHS == "" && !UNIPATH ) 
     {    cout << "You must specify PATHS unless in UNIPATH mode" << endl;
          exit(1);    }
     if (UNIPATH) 
     {    if( !COPY_NUMBER_UNKNOWN && ! IsRegularFile( copy_pdf_file ) ) 
          {    cout << "You can only CheckPredictedUnipathCoverage if the file\n  "
	            << copy_pdf_file << "\n(or specify COPY_PDF_FILE= ) has been "
                    << "generated (by UnipathCoverage)." << endl;
	       ForceAssert( IsRegularFile( copy_pdf_file ) );    }    }

     // Load unipaths (or whatever they are), and genome paths.

     if ( !PATHS_ABSOLUTE ) {
       if( UNIPATH && PATHS == "" )
	 PATHS = run_dir + "/" + READS + "." + UNIPATHS + ".k" + ToString(K);
       else 
	 PATHS = run_dir + "/" + PATHS;
     }


     // Log stream.

     String log_file = PATHS + ".PathsToLocs.out";
     
     ofstream arcout;
     if ( ARCHIVE ) arcout.open( log_file.c_str( ) );
     ostream &out = ARCHIVE? arcout : cout;
     if ( ARCHIVE ) {
       cout << "Sending output to " << log_file << "\n" << endl;
       PrintCommandPretty( out );
     }


     vecKmerPath paths(PATHS), genome_paths;
     vec<big_tagged_rpint> genome_db;
     if (RENUMBER)
     {    vecbasevector all( datadir + "/" + GENOME + ".fastb" );
          int npaths = paths.size( );
          int ngenome = all.size( );
          if ( UNIPATHS_FASTB == "" )
          {    KmerBaseBroker kbb( run_dir, K, READS );
               for ( size_t i = 0; i < paths.size( ); i++ )
                    all.push_back_reserve( kbb.Seq( paths[i] ) );    }
          else 
          {    vecbasevector unibases( run_dir + "/" + UNIPATHS_FASTB );
               all.Append(unibases);    }
          vecKmerPath universe;
          ReadsToPathsCoreY( all, K, universe, tmp_dir, NUM_THREADS );
          for ( int i = 0; i < ngenome; i++ )
               genome_paths.push_back_reserve( universe[i] );
          paths.clear( );
          for ( int i = 0; i < npaths; i++ )
               paths.push_back_reserve( universe[ngenome+i] );
          genome_db.reserve( genome_paths.sumSizes() );
          for ( size_t i = 0; i < genome_paths.size( ); i++ )
               genome_paths[i].AppendToDatabase( genome_db, i );
          Prepare(genome_db);    }
     else
     {    String genome_file = "/genome.paths.k" + ToString(K);
          if ( IsRegularFile( run_dir + genome_file ) ) 
               genome_paths.ReadAll( run_dir + genome_file );
          else if( IsRegularFile( datadir + genome_file ) ) 
          {    out << "Reading global genome.paths; I hope the kmer numbering is "
                    << "the same!" << endl;
               genome_paths.ReadAll( datadir + genome_file );    }
          else FatalErr( "No genome.paths file found!" );
          String genome_dbfile = "/genome.pathsdb_big.k" + ToString(K);
          if ( IsRegularFile( run_dir + genome_dbfile ) )
               BinaryRead2( run_dir + genome_dbfile, genome_db );
          else 
          {    genome_db.reserve( genome_paths.sumSizes() );
               for ( size_t i = 0; i < genome_paths.size( ); i++ )
                    genome_paths[i].AppendToDatabase( genome_db, i );
               Prepare(genome_db);    }    }

     // Check for canonical paths.
     
     int numNoncanonicalPaths = 0;
     for ( size_t i = 0; i < genome_paths.size(); ++i ) {
       KmerPath canonicalPath = genome_paths[i];
       canonicalPath.Canonicalize();
       if ( ! ( canonicalPath == genome_paths[i] ) )
         ++numNoncanonicalPaths;
     }
     if ( numNoncanonicalPaths > 0 ) {
       out << "Warning: " << numNoncanonicalPaths << " genome paths are not canonicalized,"
            << " so some alignments may not be found." << endl;
     }

     vec< vec<int> > genome_len( genome_paths.size( ) );
     for ( size_t i = 0; i < genome_paths.size( ); i++ )
     {    const KmerPath& p = genome_paths[i];
          genome_len[i].push_back(0);
          for ( int j = 0; j < p.NSegments( ); j++ )
               genome_len[i].push_back(
                    genome_len[i].back( ) + p.Segment(j).Length( ) );    }

     vec<genome_placement> locs;
     vec<int> actual_copy_number(paths.size(), 0);
     if( LOAD_PLACEMENTS_FROM.nonempty() ) {
       BinaryRead2( LOAD_PLACEMENTS_FROM, locs );
       for( int i=0; i < locs.isize(); i++ )
	 actual_copy_number[ locs[i].GetReadId() ] = locs[i].GetCopyNumber();
     }
     else {
       for ( size_t i = 0; i < paths.size( ); i++ )
       {  KmerPath p = paths[i];
	  int pkmers = p.KmerCount( );
          if ( pkmers < MIN_KMERS ) continue;
          if ( pkmers > MAX_KMERS ) continue;
          for ( int pass = 1; pass <= (FW_ONLY?1:2); pass++ )
          {    if ( pass == 2 ) p.Reverse( );
	       KmerPathLoc pmid = p.Begin(), pmid_copy;
	       pmid.IncrementHaltAtGap( pkmers/2 );
               longlong middle_kmer = pmid.GetKmer();
	       vec<longlong> places;
               Contains( genome_db, middle_kmer, places );
               for ( int j = 0; j < places.isize( ); j++ )
               {    const big_tagged_rpint& t = genome_db[ places[j] ];
                    int id = t.PathId( ), pp = t.PathPos( );
                    if ( id < 0 ) continue;
		    KmerPathLoc matchC(genome_paths[id],pp), matchL, matchR;
		    matchC.SetKmer(middle_kmer);
                    // These scans may fail if the genome paths are
                    // not canonicalized (which may not be possible
                    // for genome-based datasets).
                    ScanLeftPerfectMatch( pmid_copy = pmid, matchL = matchC );
                    ScanRightPerfectMatch( pmid_copy = pmid, matchR = matchC );

		    int matchlen = KmersInInterval( matchL, matchR );

		    if( (100 * matchlen)/pkmers < MATCH_LEN_PERCENT ) continue;

                    locs.push_back( genome_placement( i, pkmers, id,
			 matchL.GetLoc() + genome_len[id][matchL.GetIndex()],
			 matchR.GetLoc() + genome_len[id][matchR.GetIndex()] + K,
			 pass == 2, places.size() ) );
	       }
	  }
       }

       Sort(locs);

       for ( int i = 0; i < locs.isize( ); i++ )
	 ++actual_copy_number[ locs[i].GetReadId() ];
       for ( int i = 0; i < locs.isize( ); i++ )
	 locs[i].SetCopyNumber( actual_copy_number[ locs[i].GetReadId( ) ] );

       if( SAVE_PLACEMENTS_TO.nonempty() ) {
	 BinaryWrite2( SAVE_PLACEMENTS_TO, locs );
       }
     }

     // Compute coverage stats.

     if (STATS)
     {    int passes = ( CUTOFF == 0 ? 1 : 2 );
          for ( int pass = 1; pass <= passes; pass++ )
          {    if ( pass == 2 ) out << "With cutoff " << CUTOFF << ":\n";
               int cutoff = ( pass == 1 ? 0 : CUTOFF );
               vec< vec<Bool> > hit( genome_paths.size( ) );
               for ( size_t i = 0; i < genome_paths.size( ); i++ )
                    hit[i].resize( genome_paths[i].KmerCount( ) + K-1, False );
               for ( int i = 0; i < locs.isize( ); i++ )
               {    const genome_placement& p = locs[i];
                    if ( p.GetKmerCount( ) + K - 1 < cutoff ) continue;
                    for ( int j = p.GetStartOnGenome( ) + K/2 - 1; 
                         j <= p.GetEndOnGenome( ) - K/2; j++ )
                    {    hit[ p.GetGenomeId( ) ][j] = True;    }    }
               longlong total = 0, total_hit = 0;
               int gaps = 0;
               for( size_t i = 0; i < genome_paths.size( ); i++ )
               {    total += genome_paths[i].KmerCount( );
                    total_hit += Sum( hit[i] );
                    for ( int j = 0; j < hit[i].isize( ); j++ )
                    {    if ( !hit[i][j] ) 
                         {    ++gaps;
                              int m;
                              for ( m = j+1; m < hit[i].isize( ); m++ )
                                   if ( hit[i][m] ) break;
                              if ( pass == 1 && PRINT_GAPS )
                              {    out << "gap of size " << m-j << " at " << i 
                                        << "." << j << "-" << m
                                        << " of " << hit[i].size( ) << "\n";    }
                              j = m;    }    }    }
               vec<int> path_sizes;
               for ( size_t i = 0; i < paths.size( ); i++ )
               {    int n = paths[i].KmerCount( ) + K - 1;
                    if ( n < cutoff ) continue;
                    path_sizes.push_back(n);    }
               Sort(path_sizes);
               out << "\n" << path_sizes.size( ) << " unipaths, "
                    << "of N50 size " << N50(path_sizes) << "bp, having " << gaps
                    << " gaps, covering " << PERCENT_RATIO( 6, total_hit, total ) 
                    << " of genome.\n\n";    }

          // Compute N50 placed unipath size.

          vec<int> ulen( locs.size( ) );
          for ( int i = 0; i < locs.isize( ); i++ )
               ulen[i] = locs[i].GetKmerCount( ) + K - 1;
          Sort(ulen);
          if ( ulen.nonempty( ) )
               out << "N50 placed unipath size = " << N50(ulen) << "\n\n";    }

     // Print unplaced unipaths.

     if ( PRINT_UNPLACED < 1000000000 )
     {    vec<Bool> placed( paths.size( ), False );
          for ( int i = 0; i < locs.isize( ); i++ )
          {    const genome_placement& p = locs[i];
               placed[ p.GetReadId( ) ] = True;    }
          int count = 0;
          for ( size_t i = 0; i < paths.size( ); i++ )
               if ( !placed[i] && paths[i].KmerCount( ) >= PRINT_UNPLACED ) ++count;
          out << paths.size( ) - Sum(placed) << " unplaced unipaths\n";
          out << count << " unplaced unipaths of size >= " << PRINT_UNPLACED 
               << " kmers:\n";
          for ( size_t i = 0; i < paths.size( ); i++ )
          {    if ( !placed[i] && paths[i].KmerCount( ) >= PRINT_UNPLACED )
               {    out << "path " << i << " (" << paths[i].KmerCount( )
                         << " kmers)\n";    }    }
          out << "\n";    }

     // Load copy numbers.

     VecPdfEntryVec unipath_copy_pdf;
     if( UNIPATH ) {
       if (COPY_NUMBER_UNKNOWN) {
	 unipath_copy_pdf.resize(actual_copy_number.isize());
	 for (size_t i = 0; i < unipath_copy_pdf.size(); i++)
	   unipath_copy_pdf[i].push_back( pdf_entry( 1, 1.0 ) );
       } else
	 unipath_copy_pdf.ReadAll( copy_pdf_file );
       ForceAssertEq( actual_copy_number.size(), unipath_copy_pdf.size() );
     }

     // Show placements.

     vec<int> to_rc;
     if( SHOW_PLACEMENTS ) {

        if (RC_TOO)
        {    vec<big_tagged_rpint> pathsdb;
             CreateDatabase( paths, pathsdb );
             UnipathInvolution( paths, pathsdb, to_rc );    }
       vec<int> to_show;
       if( SHOW_THESE_CONTIGS_ONLY.nonempty() )
	 ParseIntSet( SHOW_THESE_CONTIGS_ONLY, to_show );
       int last_tig = -1, last_end_on_genome = -1;
       for ( int i = 0; i < locs.isize( ); i++ )
	 if( (locs[i].IsFw() || !SHOW_FW_ONLY)
	     &&
	     (to_show.empty() || BinMember(to_show,locs[i].GetGenomeId()))
	     &&
	     (SHOW_WINDOW_BEGIN==-1 || (SHOW_WINDOW_BEGIN <= locs[i].GetEndOnGenome() &&
					locs[i].GetStartOnGenome() <= SHOW_WINDOW_END)) ) {
           if ( last_tig == locs[i].GetGenomeId( ) 
                && last_end_on_genome - locs[i].GetStartOnGenome( ) != K-1 )
           {    out << "----------------------------------------"
                    << "--------------------------------------------\n";    }
           last_tig = locs[i].GetGenomeId( );
           last_end_on_genome = locs[i].GetEndOnGenome( );
           if ( KMER_COORDS ) {
             genome_placement& p = locs[i];
             out << "path " << setiosflags(ios::fixed) << setw(6) << p.GetReadId();
             if (RC_TOO)
             {    out << " " << setiosflags(ios::fixed) << setw(6) 
                       << to_rc[ p.GetReadId( ) ];    }
             out << resetiosflags(ios::fixed) << " (" << p.GetKmerCount() << " kmers)"
                  << " --> " << p.GetGenomeId() + TIG_OFFSET 
                  << "." << p.GetStartOnGenome() + POS_OFFSET
                  << "-" << p.GetEndOnGenome() + POS_OFFSET - K + 2 
                  << ( p.IsRc() ? " (rc)" : " (fw)" ) 
                  << " [" << p.GetCopyNumber() << " places]"
                  << "\n";
           }
           else
           {    int id = locs[i].GetReadId( );
                out << "path " << setiosflags(ios::fixed) << setw(6) << id;
                if (RC_TOO)
                {    out << " " << setiosflags(ios::fixed) << setw(6) 
                          << to_rc[ locs[i].GetReadId( ) ];    }
                int CN = -1;
                prob_t maxp = 0;
                if (UNIPATH)
                {    for ( unsigned int j = 0; j < unipath_copy_pdf[id].size( ); j++ )
                     {    if ( unipath_copy_pdf[id][j].second > maxp ) 
                          {    CN = unipath_copy_pdf[id][j].first;
                               maxp = unipath_copy_pdf[id][j].second;    }    }    }
                out << " (" << setw(4) 
                     << locs[i].GetKmerCount( ) << " kmers)" << " --> " 
                     << locs[i].GetGenomeId( ) + TIG_OFFSET 
                     << "." << locs[i].GetStartOnGenome( ) + POS_OFFSET
                     << "-" << locs[i].GetEndOnGenome( ) + POS_OFFSET
                     << ( locs[i].IsRc( ) ? " (rc)" : " (fw)" ) << " [" 
                     << locs[i].GetCopyNumber( ) << " places]" 
                     << resetiosflags(ios::fixed);
                if (UNIPATH)
                {    String prob = ToString( 100.0 * maxp, 1 );
                     if ( prob.size( ) == 4 ) prob = " " + prob;
                     out << "; CN=" << CN << "[" << prob << "%]";    }
                out << "\n";    }

         }
     }


     if( UNIPATH ) {
       // How good are the probabilities we assigned to copy numbers?
       // This requires that reads.unipaths.predicted_count.k
       // has already been generated (by UnipathCoverage).

       // Make bins based on predicted probabilities, and count how
       // often those predictions were made and were correct.

       // First process CN_TIGS argument.

       vec<Bool> excluded( paths.size( ), False );
       if ( CN_TIGS != "" )
       {    vec<int> cn_tigs;
            ParseIntSet( CN_TIGS, cn_tigs );
            for ( int i = 0; i < locs.isize( ); i++ )
            {    if ( !BinMember( cn_tigs, locs[i].GetGenomeId( ) ) )
                      excluded[ locs[i].GetReadId( ) ] = True;    }    }

       int num_bins = 20, bin, br;
       vec<longlong> numer(num_bins+1,0), denom(num_bins+1,0);

       // Break down error rates into predictions of 0, 1, 2, 3+

       vec< vec<longlong> > numer_break(4,numer);
       vec< vec<longlong> > denom_break(4,denom);

       for( size_t i=0; i < paths.size(); i++ ) {
         if ( excluded[i] ) continue;
	 int pkmers = paths[i].KmerCount();
	 if( pkmers < MIN_KMERS ) continue;
	 if( pkmers > MAX_KMERS ) continue;
	 const PdfEntryVec& pdf = unipath_copy_pdf[i];
	 int c = actual_copy_number[i];
	 double predict_max_p=0.0, predict_true_p=0.0;
	 for( unsigned int j=0; j < pdf.size(); j++ ) {
	   predict_max_p = Max( predict_max_p, pdf[j].second );
	   bin = int(floor( num_bins * pdf[j].second ));
	   br = ( pdf[j].first > 2 ? 3 : pdf[j].first );
	   denom[bin]++; denom_break[br][bin]++;
	   if( pdf[j].first == c ) {
	     numer[bin]++; numer_break[br][bin]++;
	     predict_true_p = pdf[j].second;
	   }
	 }
// 	 if( predict_max_p > 10*predict_true_p ) {
// 	   out << "Unipath " << i << ", " << pkmers << " kmers"
// 		<< ": actually " << c << ", predicted ";
// 	   for( int j=0; j < pdf.size(); j++ )
// 	     out << pdf[j].first << "(" << setprecision(3) << pdf[j].second << ") ";
// 	   out << endl;
// 	 }
       }

       out << "\n\n bin accuracy pred=0 pred=1 pred=2 pred3+  tot#predictions\n";
       for(bin = 0; bin <= num_bins; bin++) {
	 out << setw(3) << (100 * bin) / num_bins << "%   ";

	 if( denom[bin] == 0 )
	   out << "  --    ";
	 else
	   out << setw(3) << (100*numer[bin])/denom[bin] << "%    ";

	 for(int br=0; br < denom_break.isize(); br++) {
	   if( denom_break[br][bin] == 0 )
	     out << "  --   ";
	   else if( denom_break[br][bin] < 10 )
	     out << " " << numer_break[br][bin] 
		  << "/" << denom_break[br][bin] << "   ";
	   else if( denom_break[br][bin] < 100 && numer_break[br][bin]<10 )
	     out << numer_break[br][bin] 
		  << "/" << denom_break[br][bin] << "   ";
	   else
	     out << setw(3) << (100*numer_break[br][bin])/denom_break[br][bin] 
		  << "%   ";
	 }
	 out << "  " << denom[bin] << "\n";
       }
       out << endl;

       if( UNIPATH_PLOT ) {
	 using namespace ns_whiteboard;
	 whiteboard board;

	 const float x_scale = 0.01;
	 const float y_scale = 25.0;
	 const float y_offset = 10 * y_scale;
	 
	 // A baseline for each contig:
	 float x_max=0, y_max=0;

	 for( int i=0; i<genome_len.isize(); i++ ) {
	   float x1 = x_scale * genome_len[i].back();
	   board.Add(new line( xy_coords( 0, i*y_offset ),
			       xy_coords( x1, i * y_offset ),
			       2.0, black ));
	   x_max = Max(x_max, x1);
	 }

	 // For each placement of a unipath, plot its true and guessed copy number.

	 for( vec<genome_placement>::iterator pl = locs.begin(); pl != locs.end(); 
              pl++ ) {
	   int true_copy_no = pl->GetCopyNumber();
	   const PdfEntryVec& pdf = unipath_copy_pdf[pl->GetReadId()];
	   float guess_avg = 0;
	   for(unsigned int j=0; j<pdf.size(); j++)
	     guess_avg += pdf[j].first * pdf[j].second;
	   
	   float
	     x0 = x_scale * pl->GetStartOnGenome(), 
	     x1 = x_scale * pl->GetEndOnGenome(),
	     y_true = y_offset * pl->GetGenomeId() + y_scale * true_copy_no,
	     y_guess = y_offset * pl->GetGenomeId() + y_scale * guess_avg;

	   board.Add(new line( xy_coords( x0, y_guess ),
			       xy_coords( x1, y_guess ),
			       2.0, red ));
	   board.Add(new line( xy_coords( x0, y_true ),
			       xy_coords( x1, y_true ),
			       1.0, green ));
	   y_max = Max(y_max, Max(y_true, y_guess));
	 }

	 String ps_filename = "unipath_copy_number.ps";
	 Ofstream( ps_file, ps_filename );
	 ps_display ps_out( ps_file, x_max, y_max, 10.0 );
	 board.DisplayOn( &ps_out );
	 out << "\nPostscript plot saved in " << ps_filename << endl;

	 board.DeletePointers();

       } // if( UNIPATH_PLOT )

     } // if( UNIPATH )

     // Compute and print the unipath adjacency graph and predicted copy numbers.

     if (PRINT_ADJACENCY_GRAPH)
     {    String KS = ToString(K);
          vecKmerPath paths( PATHS.Before( "unipaths" ) + "paths.k" + ToString(K) ); 
          vecKmerPath paths_rc( PATHS.Before( "unipaths" )
               + "paths_rc.k" + ToString(K) ); 
          BREAD2( PATHS.Before( "unipaths" )
               + "pathsdb.k" + KS, vec<tagged_rpint>, pathsdb );
          vecKmerPath unipaths( PATHS );
          BREAD2( PATHS.Before( "unipaths" ) + "unipaths"
               + "db.k" + KS, vec<tagged_rpint>, unipathsdb );
          digraph A;
          BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths,
               unipathsdb, A );
          out << "\nAdjacency graph and predicted copy numbers:\n";
          for ( int v = 0; v < A.N( ); v++ )
          {    out << "\n" << v;
               if ( SHOW_PLACEMENTS && RC_TOO ) out << "[" << to_rc[v] << "]";
               if ( !COPY_NUMBER_UNKNOWN )
               {
               out << ":";
               for ( unsigned int j = 0; j < unipath_copy_pdf[v].size( ); j++ )
               {    out << " " << unipath_copy_pdf[v][j].first << "(" << setprecision(3)
                         << unipath_copy_pdf[v][j].second << ")";    }
               }
               out << "\n";
               out << v << " <--";
               for ( int j = 0; j < A.To(v).isize( ); j++ )
                    out << " " << A.To(v)[j];
               out << "\n";
               out << v << " -->";
               for ( int j = 0; j < A.From(v).isize( ); j++ )
                    out << " " << A.From(v)[j];
               out << "\n";    }    }

}
