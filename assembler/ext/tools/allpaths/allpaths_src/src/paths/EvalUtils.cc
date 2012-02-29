///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "Equiv.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/AlignHyperKmerPath.h"
#include "paths/EvalUtils.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/ReadsToPathsCoreX.h"

// MakeDepend: dependency ScaffoldAccuracy

/**
   Function: SummarizeReferenceCoverage

   Print out information reporting gaps in coverage (tiny gaps are ignored).
   
*/
void SummarizeReferenceCoverage
( longlong& total_bases, longlong& total_covered, ostream& out,
  const vecbasevector& genome, const vecbasevector& genome_diploid,
  const vecbitvector& genome_amb, const vec<look_align>& aligns,
  const Bool brief )
{    size_t n_ref_contigs = genome.size( );
     vec< vec<ho_interval> > covered( n_ref_contigs );
     longlong ambig = 0;
     if ( genome_amb.size( ) > 0 )
     {    for ( size_t i = 0; i < n_ref_contigs; i++ )
          {    for ( unsigned int j = 0; j < genome[i].size( ); j++ )
               {    if ( !genome_amb[i][j] ) continue;
                    int k = genome_amb[i].NextDiff(j);
                    covered[i].push( j, k );
                    ambig += k - j;
                    j = k - 1;    }    }    }
     for ( int i = 0; i < aligns.isize( ); i++ )
     {    const look_align& la = aligns[i];
          covered[la.target_id].push( la.pos2( ), la.Pos2( ) );    }
     total_bases = 0, total_covered = 0;
     vec<int> gaps;
     for ( size_t i = 0; i < n_ref_contigs; i++ )
     {    vec< pair<ho_interval, int> > covered2;
          CondenseIntervals( genome[i].size( ), covered[i], covered2 );
          for ( unsigned int j = 0; j < covered2.size( ); j++ )
          {    if ( covered2[j].second == 0 )
                    gaps.push_back( covered2[j].first.Length( ) );    }    }
     ReverseSort(gaps);
     int min_gap = 1, gaps_to_show = 25;
     if ( gaps.isize( ) > gaps_to_show ) min_gap = gaps[gaps_to_show-1];
     if ( !brief )  out << "\ngaps of size >= " << min_gap << ":\n";
     for ( size_t i = 0; i < n_ref_contigs; i++ )
     {    vec< pair<ho_interval, int> > covered2;
          CondenseIntervals( genome[i].size( ), covered[i], covered2 );
          for ( unsigned int j = 0; j < covered2.size( ); j++ )
          {    total_bases += covered2[j].first.Length( );
               if ( covered2[j].second > 0 )
                    total_covered += covered2[j].first.Length( );
               else
               {    if ( !brief && covered2[j].first.Length( ) >= min_gap )
                    {    out << i << "." << covered2[j].first.Start( ) << "-" 
                              << covered2[j].first.Stop( ) << " (" 
                              << covered2[j].first.Length( )
                              << " bases)\n";    }    }    }    }
     total_covered -= ambig;
     total_bases -= ambig;
     out << "\nSUMMARY: " << PERCENT_RATIO(4, total_covered, total_bases)
	 << " (" << total_covered << "/" << total_bases << ") "
	 << "of reference covered (excluding ambiguous bases)" << endl;
     
     // If a diploid version of the genome (from the file genome.diploid.fastb,
     // as used in SimulateReads) was supplied, use it to perform some further
     // analysis.
     if ( genome_diploid.size( ) > 0 ) {
       ForceAssertEq( genome_diploid.size( ), 2 * n_ref_contigs );
       
       // We assume the diploid genome is identical to the original genome,
       // doubled, with possible mutations but no indels.  This is a safe
       // assumption for genome.diploid.fastb files created by MutateReference.
       for ( size_t i = 0; i < n_ref_contigs; i++ ) {
	 ForceAssert( genome[i] == genome_diploid[i] );
	 ForceAssert( genome_diploid[i].size( ) ==
		      genome_diploid[i + n_ref_contigs].size( ) );
       }
       
       // Find all SNPs and remove those loci from the coverage statistics.
       longlong SNPs = 0, SNPs_covered = 0;
       for ( size_t i = 0; i < n_ref_contigs; i++ )
	 for ( unsigned int j = 0; j < genome[i].size( ); j++ )
	   if ( genome_diploid[i][j] != genome_diploid[i + n_ref_contigs][j] ) {
	     
	     // We have found a SNP!
	     SNPs++;
	     // Determine if this SNP is covered already.
	     // This is a dumb-ass algorithm that might need to be optimized;
	     // it's O(N^2) when it could be O(N log N).
	     for ( int k = 0; k < covered[i].isize( ); k++ )
	       if ( covered[i][k].Contains( j ) ) {
		 SNPs_covered++;
		 break;
	       }
	   }
       
       // Report.
       total_covered -= SNPs_covered;
       total_bases -= SNPs;
       out << "\nSUMMARY: " << PERCENT_RATIO(4, total_covered, total_bases)
	   << " (" << total_covered << "/" << total_bases << ") "
	   << "of reference covered (excluding ambiguous bases and SNPs)" << endl;
     }
}

/**
   Function: ComputeComponentSizes
   
   Compute the component sizes.  The size of each
   component here is defined to be the sum of the component's
   edges, which is not quite right but usually close enough.  To
   compensate for duplicate edges, we only take the longest edge
   between any two vertices.

   Output params:

     component_sizes - the size of each component
   
 */
void ComputeComponentSizes( const HyperKmerPath& h, vec< int >& component_sizes ) {


  equiv_rel e;
  h.ComponentRelation(e);
  vec<int> reps;
  e.OrbitRepsAlt(reps);
  component_sizes.resize( reps.size( ), 0 );
  for ( int i = 0; i < reps.isize( ); i++ )
    {
      vec<int> o;
      e.Orbit( reps[i], o );
      for ( int j = 0; j < o.isize( ); j++ )
	{    int v = o[j];
	  map<int,int> longestEdgeByVertex;
	  for ( int l = 0; l < h.From(v).isize( ); l++ ) {    
	    int vertex = h.From(v)[l];
	    int length = h.EdgeObjectByIndexFrom( v, l ).KmerCount();
	    map<int,int>::iterator exists = longestEdgeByVertex.find( vertex );
	    if ( exists != longestEdgeByVertex.end() )
	      exists->second = max( exists->second, length );
	    else
	      longestEdgeByVertex.insert( make_pair( vertex, length ) );
	  }
	  for ( map<int,int>::iterator longest = longestEdgeByVertex.begin();
		longest != longestEdgeByVertex.end(); ++longest )
	    component_sizes[i] += longest->second;
	}
    }
}  // ComputeComponentSizes()

/**
    Function: PrintGraphStatistics

    Print out summary stats about the HyperKmerPath representing the entire
    assembly.

*/
void PrintGraphStatistics
( ostream& out, const HyperKmerPath& h, 
  const vecbasevector& genome, const vecbasevector& genome_diploid,
  const vec<look_align>& aligns, const vec< vec<int> >& aligns_index )
{
     out << "\nGRAPH STATISTICS:\n";
     int sourcesink = 0, ncomp = h.ConnectedComponents( );
     for ( int v = 0; v < h.N( ); v++ )
     {    if ( h.Source(v) ) ++sourcesink;
          if ( h.Sink(v) ) ++sourcesink;    }    

     vec< int > component_sizes;
     ComputeComponentSizes( h, component_sizes );
     ForceAssertEq( ncomp, component_sizes.isize() );
     Sort( component_sizes );
     longlong csum = BigSum( component_sizes);
     
     out << ncomp << " components";
     if ( ncomp > 0 ) out << ", of N50 size " << N50(component_sizes) << " (total size " << csum << ")";
     out << "\n";
     size_t n_ref_contigs = genome.size( );
     if ( n_ref_contigs != 0 ) { // i.e., if a reference genome was provided
       vec<unsigned int> reference_sizes;
       reference_sizes.reserve( n_ref_contigs );
       for ( size_t i = 0; i < n_ref_contigs; ++i )
	 reference_sizes.push_back( genome[i].size() );
       Sort( reference_sizes );
       out << "(reference has " << n_ref_contigs << " sequence"
	   << ( n_ref_contigs > 1 ? "s, with N50 size " : ", of size " )
	   << N50(reference_sizes) << ")\n";
     }

     // Describe vertices and edges.

     out << h.N( ) << " vertices\n";
     int nedges = h.EdgeObjectCount( );
     out << nedges << " edges, of N50 size " << h.EdgeN50( ) << "\n";

     // Describe ambiguities.

     int amb = ncomp + nedges - h.N( ), loop1 = 0, loop2 = 0, loop34 = 0;
     for ( int v = 0; v < h.N( ); v++ )
     {    for ( int j = 0; j < h.From(v).isize( ); j++ )
          {    int w = h.From(v)[j];
               if ( w != v ) continue;
               const KmerPath& p = h.EdgeObjectByIndexFrom( v, j );
               int n = p.KmerCount( );
               if ( n == 1 ) ++loop1;
               else if ( n == 2 ) ++loop2;    
               else if ( n == 3 || n == 4 ) ++loop34;    }    }
     
     // If a diploid version of the genome (from the file genome.diploid.fastb,
     // as used in SimulateReads) was supplied, we can classify some of the
     // ambiguities as diploid, i.e., arising from SNPs.
     int n_dips = 0;
     if ( genome_diploid.size( ) > 0 ) {
       ForceAssertEq( genome_diploid.size( ), 2 * n_ref_contigs );
       
       // We assume the diploid genome is identical to the original genome,
       // doubled, with possible mutations but no indels.  This is a safe
       // assumption for genome.diploid.fastb files created by MutateReference.
       for ( size_t i = 0; i < n_ref_contigs; i++ ) {
	 ForceAssert( genome[i] == genome_diploid[i] );
	 ForceAssert( genome_diploid[i].size( ) ==
		      genome_diploid[i + n_ref_contigs].size( ) );
       }
       
       vec<int> SNP_vertices;
       
       for ( int v = 0; v < h.N( ); v++ ) {
	 // Look for vertices v that are located at the opening bit of a bubble:
	 //        _-_
	 //  --> v     w -->
	 //        -_-
	 // (there may be more than one edge entering v or leaving w)
	 if ( h.FromSize(v) != 2 ) continue;
	 if ( h.From(v)[0] != h.From(v)[1] ) continue;
	 
	 // Find the alignments of the two edges in the bubble.
	 // Require them to align uniquely to the reference, if at all
	 // (the SNPy edge might not have aligned.)
	 int e1 = h.EdgeObjectIndexByIndexFrom( v, 0 );
	 int e2 = h.EdgeObjectIndexByIndexFrom( v, 1 );
	 if ( aligns_index[e1].size( ) > 1 ) continue;
	 if ( aligns_index[e2].size( ) > 1 ) continue;
	 if ( aligns_index[e1].empty( ) && aligns_index[e2].empty( ) ) continue;
	 const look_align & la = aligns_index[e1].size( ) > 0 ?
	   aligns[ aligns_index[e1][0] ] :
	   aligns[ aligns_index[e2][0] ];
	 
	 // Look for a SNP on the region of genome covered by this alignment.
	 // Ignore the first and the last K-1 bases of the covered region;
	 // no need to deal with the alignment's orientation.
	 int K = h.K( );
	 int genome_ID = la.TargetId( );
	 int genome_start_loc = la.StartOnTarget( ) + K - 1;
	 int genome_end_loc   = la.EndOnTarget( )   - K + 1;
	 
	 for ( int i = genome_start_loc; i < genome_end_loc; i++ ) {
	   // If we find a SNP here, we mark this as a diploid ambiguity.
	   if ( genome_diploid[genome_ID][i] != genome_diploid[genome_ID + n_ref_contigs][i] ) {
	     n_dips++;
	     SNP_vertices.push_back( v );
	     break;
	   }
	 }
       }
       
       // Now that we've found all diploid ambiguities, we create a version
       // of the HyperKmerPath with those ambiguities removed, just for the
       // sake of finding its graph statistics.
       HyperKmerPath hkp_dip = h;
       hkp_dip.PopBubbles( SNP_vertices );
       
       out << "(Ignoring diploidy, graph has " << hkp_dip.N()
	   << " vertices, and " << hkp_dip.EdgeObjectCount()
	   << " edges with N50 size " << hkp_dip.EdgeN50() << ".)" << endl;
     }
       
     out << "ambiguities:\n";
     if (genome_diploid.size( ) > 0) out << "    " << n_dips << " diploid ambiguities\n";
     out << "    " << loop1 << " perfect mononukes\n";
     out << "    " << loop2 << " perfect dinukes\n";
     out << "    " << loop34 << " perfect trinukes or tetranukes\n";
     out << "    " << amb - loop1 - loop2 - loop34 - n_dips << " other\n";    
     out << "    " << "-------------------------------------\n";
     out << "    " << amb - n_dips << " total";
     if (genome_diploid.size( ) > 0) out << " (non-diploid)";
     out << "\n";
     if ( n_ref_contigs != 0 ) {
       double genome_size = genome.sumSizes();
       double amb_per_Mb = 
         static_cast<double>(amb - n_dips)/static_cast<double>(genome_size)*1000000.;
       out << ToString( amb_per_Mb ) << " ambiguities per Mb of reference\n";
     }
}

/**
   Function: EvaluateAssembly

   Print out an evaluation of the assembly, using truth data if
   available.  Note that this will reorder the HyperKmerPath h.

   Called at the end of <LocalizeReadsTail()>.

   Input/output parameters:

      h - a <HyperKmerPath> representing *the entire assembly*.
        On output, its <components> are reordered to follow the
	reference.

   NOTE: THIS FUNCTION IS DEPRECATED.  It is used only in the original
   (non-LG) RunAllPaths pipeline, and it is not maintained for LG use.

*/
void EvaluateAssembly( HyperKmerPath& h, const KmerBaseBroker* kbb,
     const String& data_dir, const String& wrun_dir, const String& sub_dir,
     const vecbasevector& genome, const vecbitvector& genome_amb, 
     const Bool DIPLOID, const Bool USE_TRUTH, const Bool FILTER_ALIGNS, 
     const Bool WRITE_ALIGNS, const Bool REORDER, nbases_t MIN_TRUSTED_PATH_LEN,
     const filename_t& HYPER, String report_suffix )
{   
     if ( report_suffix != "" ) report_suffix = "." + report_suffix;

     // Align merged HyperKmerPath edges to reference, then reorder the 
     // components so that they follow the reference.

     vec<look_align> aligns;
     vec< vec<int> > aligns_index;

     vecbasevector genome_diploid;
     
     vec<TrustedPath> trusted_paths;
     if (USE_TRUTH)
     {    AlignHyperKmerPath( h, kbb, data_dir + "/genome", wrun_dir, aligns,
               aligns_index );

          if (FILTER_ALIGNS) 
               FilterAligns( h, aligns, aligns_index, trusted_paths, MIN_TRUSTED_PATH_LEN );


	  if (REORDER) {
            double reorderclock = WallClockTime();
            ReorderToFollowReference( h, aligns, aligns_index );
            cout << TimeSince(reorderclock) 
                 << " used reordering to follow reference." << endl;
            if (FILTER_ALIGNS) 
              FilterAligns( h, aligns, aligns_index, trusted_paths, MIN_TRUSTED_PATH_LEN );
	  }

          if (WRITE_ALIGNS) 
          {    double writeclock = WallClockTime();
               WriteLookAlignBinary( sub_dir + "/" + HYPER + ".aligns", aligns ); 
               cout << TimeSince(writeclock) 
                    << " used writing aligns." << endl;    }

          BinaryWrite( sub_dir + "/" + HYPER + ".trusted_paths", trusted_paths );
     }

     int K = h.K( );
     vec<int> to_right_vertex, to_left_vertex;
     h.ToRight(to_right_vertex), h.ToLeft(to_left_vertex);

     // Print text version of merged HyperKmerPath.

     String bar = "=============================================================="
          "======================";
     String dbar = ".............................................................."
          "......................";
     if (USE_TRUTH)
     {    cout << "\n" << bar << "\n"
               << "EVALUATION OF ASSEMBLY VERSUS REFERENCE\n" << bar << "\n";
          PrintAlignedHyperKmerPath( cout, h, kbb, genome, aligns, 
               aligns_index, True, &trusted_paths );
          if (DIPLOID) // bad, but not getting done elsewhere
          {    Ofstream( bout, sub_dir + "/report.brief" + report_suffix );
               bout << "EVALUATION OF ASSEMBLY VERSUS REFERENCE\n";
               PrintAlignedHyperKmerPath( bout, h, kbb, genome, aligns, 
                    aligns_index, True, &trusted_paths );    
               ReportMisassemblies( bout, h, aligns, aligns_index );
               longlong total_bases, total_covered;
               SummarizeReferenceCoverage( total_bases, total_covered, bout, genome, genome_diploid, genome_amb, aligns, False );
               int mismatches = 0, indels = 0;
               for ( int i = 0; i < aligns.isize( ); i++ )
               {    mismatches += aligns[i].mutations;
                    indels += aligns[i].indels;    }
               bout << "\nBASE ERRORS: " << mismatches << " mismatches and "
                    << indels << " indels.\n";
               PrintGraphStatistics( bout, h, genome, genome_diploid, aligns, aligns_index );    
               if ( IsRegularFile( sub_dir + "/linear_scaffolds.fasta" ) )
               {    bout << "\nSCAFFOLD STATISTICS (excluding contigs < 1 kb):\n";
                    vec<String> report = AllOfOutput( "ScaffoldAccuracy" 
                         + ARG(NH, True) + ARG(REFHEAD, data_dir + "/../genome")
                         + ARG(MIN_CONTIG, 0)
                         + ARG(ASSEMBLY, sub_dir + "/linear_scaffolds.fasta") );
                    for ( int i = 0; i < report.isize( ); i++ )
                         bout << report[i] << "\n";    }    }    }
     else h.PrintSummaryPlus( cout, 0, 0, 0, 0, 0, False );

     // Write brief and full reports in haploid case.

     if ( USE_TRUTH && !DIPLOID )
     {    for ( int pass = 1; pass <= 2; pass++ )
          {    Ofstream( bout, sub_dir + "/report." 
                    + ( pass == 1 ? "brief" : "full" ) + report_suffix );
               bout << "EVALUATION OF ASSEMBLY VERSUS REFERENCE\n";
               PrintAlignedHyperKmerPath( bout, h, kbb, genome, aligns, 
                    aligns_index, True, &trusted_paths, pass == 1, DIPLOID );
               bout << "\n" << dbar << "\n";
               longlong total_bases, total_covered;
               vecbitvector genome_amb;
               SummarizeReferenceCoverage( total_bases, total_covered, bout, 
					   genome, genome_diploid, genome_amb, aligns, pass == 1 );
               if ( pass == 2 )
               {    bout << "\nEdges of size >= 10 kb that do not align "
                         << "to reference at all:\n\n";
                    vec<Bool> aligned( h.EdgeObjectCount( ), False );
                    for ( int i = 0; i < aligns.isize( ); i++ )
                         aligned[ aligns[i].query_id ] = True;
                    longlong total_covered_plus = total_covered;
                    for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
                    {    int n = h.EdgeObject(i).KmerCount( );
                         if ( aligned[i] ) continue;
                         total_covered_plus += n;
                         if ( n < 10000 ) continue;
                         bout << "[" << BaseAlpha(i) << "] - " 
                              << n << " kmers\n";    }
                    bout << "\nSUMMARY: " 
                         << PERCENT_RATIO(4, total_covered_plus, total_bases)
                         << " (" << total_covered_plus << "/" << total_bases << ") "
                         << "of reference covered (including unaligned edges)\n";
                    bout << "\nEdges of size >= 10 kb that align imperfectly "
                         << "to reference:\n\n";
                    vec<Bool> perf( h.EdgeObjectCount( ), False );
                    for ( int i = 0; i < aligns.isize( ); i++ )
                    {    if ( aligns[i].FullLength( ) && aligns[i].Errors( ) == 0 )
                              perf[ aligns[i].query_id ] = True;    }
                    for ( int i = 0; i < perf.isize( ); i++ )
                    {    int n = h.EdgeObject(i).KmerCount( );
                         if ( !aligned[i] || perf[i] || n < 10000 ) continue;
                         bout << "[" << BaseAlpha(i) << "] - " 
                              << n << " kmers\n";    }
                    vec<Bool> short_or_bad( aligns.size(), False );
                    for ( int i = 0; i < aligns.isize( ); i++ )
                      if ( aligns[i].query_length < 10000 ||
                           aligns[i].Errors() != 0 ||
                           ! aligns[i].FullLength() )
                        short_or_bad[i] = True;
                    vec<look_align> good_and_long_aligns = aligns;
                    EraseIf( good_and_long_aligns, short_or_bad );
                    bout << "\nCoverage of reference by edges of size >= 10 kb "
                         << "that align perfectly:";
                    longlong total_bases, total_covered;
                    vecbitvector genome_amb;
                    SummarizeReferenceCoverage( total_bases, total_covered, bout,
						genome, genome_diploid, genome_amb, good_and_long_aligns, True );    }
               PrintGraphStatistics( bout, h, genome, genome_diploid, aligns, aligns_index );    }    }

     // Check for putative misassemblies.

     if (USE_TRUTH) ReportMisassemblies( cout, h, aligns, aligns_index );

     // Summarize coverage of reference.

     if (USE_TRUTH) 
     {    longlong total_bases, total_covered;
          vecbitvector genome_amb;
          SummarizeReferenceCoverage( total_bases, total_covered, cout, genome, 
				      genome_diploid, genome_amb, aligns, False );
     }

     // Output graph statistics.

     PrintGraphStatistics( cout, h, genome, genome_diploid, aligns, aligns_index );

     // For diploid case, attempt to do some repairs.

     if ( DIPLOID && USE_TRUTH && genome.size( ) % 2 == 0 )
     {    cout << "\n" << bar << "\n"
               << "USING REFERENCE TO IMPROVE ASSEMBLY\n" << bar << "\n\n";
          HyperBasevector hb( h, *kbb );
          Ofstream( bout, sub_dir + "/report.brief" + report_suffix );
          bout << "EVALUATION OF MODIFIED DIPLOID ASSEMBLY VERSUS REFERENCE\n";

          // Find cases where only one allele is covered.

          vec< vec<ho_interval> > covered( genome.size( ) );
          vec< vec<int> > covered_by( genome.size( ) );
          for ( int i = 0; i < aligns.isize( ); i++ )
          {    const look_align& la = aligns[i];
               covered[la.target_id].push_back(
                    ho_interval( la.pos2( ), la.Pos2( ) ) );
               covered_by[la.target_id].push_back(i);    }

          // Now covered tells us which regions are covered by an alignment, an
          // covered_by tells us the corresponding alignment that the coverage came
          // from.

          vec<int> L( genome.size( ), 0 );
          for ( size_t g = 0; g < genome.size( ); g++ )
          {    SortSync( covered[g], covered_by[g] );
               for ( int i = 0; i < covered[g].isize( ); i++ )
                    L[g] = Max( L[g], covered[g][i].Length( ) );    }

          // Now L provides upper bounds for the lengths of the intervals in
          // covered.

          // TODO: potentially dangerous truncation of index by pair<int,int>
          vec< pair<int,int> > snps_added;
          for ( size_t g = 0; g < genome.size( )/2; g++ )
          {    vec< vec<ho_interval> > cov(2);
               vec< vec<int> > cov_by(2);
               cov[0] = covered[g];
               cov[1] = covered[ g + genome.size( )/2 ];
               cov_by[0] = covered_by[g];
               cov_by[1] = covered_by[ g + genome.size( )/2 ];
               for ( int j1 = 0; j1 < 2; j1++ )
               {    int j2 = 1 - j1;
                    size_t g1 = g, g2 = g + genome.size( )/2;
                    if ( j1 == 1 ) swap( g1, g2 );

                    // Now we have homologous chromosomes g1 and g2.

                    for ( int i1 = 0; i1 < cov_by[j1].isize( ); i1++ )
                    {    const look_align& la1 = aligns[ cov_by[j1][i1] ];
                         if ( la1.Errors( ) > 0 ) continue;
                         if ( !aligns_index[ la1.query_id ].solo( ) )
                              continue;

                         // Now we have a perfect alignment la1 to chromosome 1.
                         // Is there a SNP under it which is uncovered on 
                         // chromosome 2?

                         ho_interval h1( la1.pos2( ), la1.Pos2( ) );
                         vec<int> I;
                         OverlapIndices( h1, cov[j2], L[g2], I );
                         vec<int> snps_missing2;
                         for ( int u = la1.pos2( ); u < la1.Pos2( ); u++ )
                         {    if ( genome[g1][u] == genome[g2][u] ) continue;

                              // Now have a SNP at position u lying under la1.
                              // See if it is under an alignment on the other 
                              // chromosome.

                              Bool cov2 = False;
                              for ( int k = 0; k < I.isize( ); k++ )
                              {    const look_align& la2 
                                        = aligns[ cov_by[j2][ I[k] ] ];

                                   // Now we have an alignment la2 to chromosome 2,
                                   // whose coverage overlaps that of la1.
                                   //
                                   // See if u lies under la2.  Check for the 
                                   // special case where u lies in the K-1 end bases 
                                   // and has the wrong base.

                                   if ( !Member( la2.Extent2( ), u ) ) continue;
                                   if ( ( u < la1.pos2( ) + (K-1)
                                          || u >= la1.Pos2( ) - (K-1) )
                                        && 
                                        ( u < la2.pos2( ) + (K-1) 
                                          || u >= la2.Pos2( ) - (K-1) ) 
                                        && 
                                        !la2.rc1 )
                                   {    const basevector& e2 = 
                                             hb.EdgeObject( la2.query_id );
                                        if ( genome[g2][u] != e2[ u - la2.pos2( ) ] )
                                             continue;    }

                                   // Now we know that u lies under la2, so we don't
                                   // want to edit the assembly.

                                   cov2 = True;
                                   break;    }
                              if (cov2) continue;

                              // Now we know that u is covered by the perfect 
                              // alignment la1 on chromosome 1, but is not covered
                              // at all on chromosome 2.
          
                              cout << "see missing SNP at base " << g2 << "."
                                   << u << endl;
                              snps_added.push( g2, u );
                              snps_missing2.push_back(u);    }
                         if ( snps_missing2.empty( ) ) continue;

                         // Generate new edge.

                         basevector e;
                         e.SetToSubOf( genome[g2], la1.pos2( ), la1.extent2( ) );
                         int id = la1.query_id;
                         int v = to_left_vertex[id], w = to_right_vertex[id];
                         cout << "adding edge from " << v << " to " << w
                              << " having " << snps_missing2.size( ) << " SNPs\n";
                         hb.AddEdge( v, w, e );    }    }    }
          int total_snps = 0;
          for ( size_t g = 0; g < genome.size( )/2; g++ )
          {    for ( unsigned int j = 0; j < genome[g].size( ); j++ )
                    if ( genome[g][j] != genome[ g + genome.size( )/2 ][j] )
                         ++total_snps;    }
          UniqueSort(snps_added);
          cout << "\n" << snps_added.size( ) << " SNPS ADDED ("
               << PERCENT_RATIO( 4, snps_added.isize( ), total_snps ) 
               << ")\n" << endl;

          // Build HyperKmerPath h2 corresponding to hb.  This is really ugly
          // and stupid.  We ought to instead make the subsequent alignment, etc.
          // take a HyperBasevector hb as input.

          vecbasevector bases;
          for ( int j = 0; j < hb.EdgeObjectCount( ); j++ )
               bases.push_back( hb.EdgeObject(j) );
          Mkdir777( wrun_dir + "/2" );
          bases.WriteAll( wrun_dir + "/2/reads.fastb" );
          vecKmerPath spaths;
          String KS = ToString(K);
          ReadsToPathsCoreY( bases, K, spaths );
          spaths.WriteAll( wrun_dir + "/2/reads.paths.k" + KS );
          vecKmerPath spaths_rc(spaths);
          for ( size_t i = 0; i < spaths.size( ); i++ )
               spaths_rc[i].Reverse( );
          spaths_rc.WriteAll( wrun_dir + "/2/reads.paths_rc.k" + KS );
          vec<tagged_rpint> spathsdb;
          CreateDatabase( spaths, spaths_rc, spathsdb );
          BinaryWrite2( wrun_dir + "/2/reads.pathsdb.k" + KS, spathsdb );
          KmerBaseBroker* kbb2 = new KmerBaseBroker( wrun_dir + "/2", K );
          vec<KmerPath> these_paths;
          for ( int j = 0; j < hb.EdgeObjectCount( ); j++ )
               these_paths.push_back( spaths[j] );
          HyperKmerPath h2( K, hb, these_paths );

          // Redo evaluation.

          AlignHyperKmerPath( h2, kbb2, data_dir + "/genome", wrun_dir, aligns,
               aligns_index );
          if (FILTER_ALIGNS) FilterAligns( h2, aligns, aligns_index, trusted_paths, MIN_TRUSTED_PATH_LEN );
          ReorderToFollowReference( h2, aligns, aligns_index );
          cout << "\n" << bar << "\n"
               << "EVALUATION OF IMPROVED ASSEMBLY VERSUS REFERENCE\n" 
               << bar << "\n";
          PrintAlignedHyperKmerPath( cout, h2, kbb2, genome, aligns, 
               aligns_index, True, &trusted_paths );
          PrintAlignedHyperKmerPath( bout, h2, kbb2, genome, aligns, 
               aligns_index, True, &trusted_paths, True, DIPLOID );
          bout << "\n" << dbar << "\n\n" << snps_added.size( ) << " SNPS ADDED ("
               << PERCENT_RATIO( 4, snps_added.isize( ), total_snps ) << ")\n";
          longlong total_bases, total_covered;
          vecbitvector genome_amb;
          SummarizeReferenceCoverage( total_bases, total_covered, cout, genome, 
				      genome_diploid, genome_amb, aligns, False );
          SummarizeReferenceCoverage( total_bases, total_covered, bout, genome, 
				      genome_diploid, genome_amb, aligns, True );
          vec<Bool> short_or_bad( aligns.size(), False );
          for ( int i = 0; i < aligns.isize( ); i++ )
            if ( aligns[i].query_length < 10000 ||
                 aligns[i].Errors() != 0 ||
                 ! aligns[i].FullLength() )
              short_or_bad[i] = True;
          vec<look_align> good_and_long_aligns = aligns;
          EraseIf( good_and_long_aligns, short_or_bad );
          bout << "\nCoverage of reference by edges of size >= 10 kb that "
               "align perfectly:";
          SummarizeReferenceCoverage( total_bases, total_covered, bout, genome, 
				      genome_diploid, genome_amb, good_and_long_aligns, True );
          PrintGraphStatistics( cout, h2, genome, genome_diploid, aligns, aligns_index );
          PrintGraphStatistics( bout, h2, genome, genome_diploid, aligns, aligns_index );
          cout << "\nFixed dot file is in\n" << RealPath(sub_dir)
               << "/hyper.fixed.dot.\n";
          Ofstream( dot, sub_dir + "/hyper.fixed.dot" );
          h2.PrintSummaryDOT0w(dot);    

	  if ( IsRegularFile( sub_dir + "/linear_scaffolds.fasta" ) )
	    {    bout << "\nSCAFFOLD STATISTICS (excluding contigs < 1 kb):\n";
	    vec<String> report = AllOfOutput( "ScaffoldAccuracy" 
	         + ARG(NH, True) + ARG(REFHEAD, data_dir + "/../genome")
                 + ARG(MIN_CONTIG, 0)
		 + ARG(ASSEMBLY, sub_dir + "/linear_scaffolds.fasta") );
	    for ( int i = 0; i < report.isize( ); i++ )
	      bout << report[i] << "\n";    }     }
     if ( USE_TRUTH ) cout << "\n" << "Brief report is in\n" << RealPath(sub_dir) 
          << "/report.brief" + report_suffix + ".\n";    }

void FilterAligns( const HyperKmerPath& h, vec<look_align>& aligns,
     vec< vec<int> >& aligns_index, vec<TrustedPath>& trusted_paths,
     const int MIN_TRUSTED_PATH_LEN )
{    //double refclock = WallClockTime( );
     FilterByReference( h, h.K( ), aligns, aligns_index, trusted_paths );
     //cout << TimeSince(refclock) << " used creating trusted paths" << endl;

     if ( MIN_TRUSTED_PATH_LEN > 0 ) {
       //double lenclock = WallClockTime( );
       FilterPathsByLength( trusted_paths, MIN_TRUSTED_PATH_LEN, 0 );
       //cout << TimeSince(lenclock) << " used filtering trusted paths by length" << endl;
     }

     //double alignclock = WallClockTime( );
     FilterPathsByAlignDominance( trusted_paths );
     //cout << TimeSince(alignclock) << " used filtering trusted paths by align dominance." << endl;

     //double edgeclock = WallClockTime( );
     FilterPathsByEdgeDominance( trusted_paths, aligns_index.size() );
     //cout << TimeSince(edgeclock) << " used filtering trusted paths by edge dominance." << endl;

     TrustedPathsToIndexedAligns( trusted_paths, h.EdgeObjectCount( ),
          aligns, aligns_index );    }



