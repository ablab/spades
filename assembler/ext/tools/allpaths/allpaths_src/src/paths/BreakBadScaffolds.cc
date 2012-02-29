///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// BreakBadScaffolds.  

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "paths/BreakBadScaffolds.h"
#include "paths/ScaffoldsUtils.h"

// MakeDepend: dependency Fastb

void break_scaffolds( vec<superb> & scaffolds, vec<superb> & new_scaffolds,
		      PairsManager & pairs, vec<alignlet> & aligns0, vec<int> & aligns0_index, 
		      const int nreads, const int min_reach_away, vec<int> & trace_ids,
		      ofstream & log, vec<fastavector> & contigs, ofstream & sout, ofstream & cgout ){
 
  // Define heuristic constants.
  
  const int min_dist_from_end = 1000;
  const double dev_mult = 4.0;
  const int blink_dev_mult = 1;
  const int min_scaffold = 1000;
  const int trim_back = 80;
  const int min_contig = 1000;
  const int min_spread = 10;


  int ntigs = 0;
  for (int ii=0; ii<scaffolds.isize( ); ii++)
    for (int jj=0; jj<scaffolds[ii].Ntigs( ); jj++)
      ntigs = Max( ntigs, scaffolds[ii].Tig( jj ) );
  ntigs++;
  
  vec<int> tig_sizes( ntigs, 0 );
  for (int ii=0; ii<scaffolds.isize( ); ii++)
    for (int jj=0; jj<scaffolds[ii].Ntigs( ); jj++)
      tig_sizes[ scaffolds[ii].Tig( jj ) ] = scaffolds[ii].Len( jj );
   
   // Index scaffolds.
     cout << Date( ) << ": indexing scaffolds" << endl;
     vec<int> to_scaffold(ntigs, -1), to_scaffold_pos(ntigs, -1);
     for ( int i = 0; i < scaffolds.isize( ); i++ ) {
       const superb& s = scaffolds[i];
       for ( int j = 0; j < s.Ntigs( ); j++ ) {
	 ForceAssertEq( to_scaffold[ s.Tig(j) ], -1 );
	 to_scaffold[ s.Tig(j) ] = i;
	 to_scaffold_pos[ s.Tig(j) ] = j;
       }
     }

     // Find pairs that land logically on a scaffold, and mark the physically
     // covered parts.

     cout << Date( ) << ": finding coverage" << endl;
     size_t nbatches = 10 * omp_get_max_threads( );
     vec< vec< vec<ho_interval> > > covi(nbatches);
     for ( size_t i = 0; i < nbatches; i++ )
          covi[i].resize(ntigs);
     #pragma omp parallel for
     for ( size_t bi = 0; bi < nbatches; bi++ )
     {    size_t start = ( bi * nreads ) / nbatches;
          size_t stop = ( (bi+1) * nreads ) / nbatches;
          for ( size_t id1 = start; id1 < stop; id1++ )
          {    if ( !pairs.isPaired(id1) ) continue;
               int id2 = pairs.getPartnerID(id1);
               if ( aligns0_index[id1] < 0 || aligns0_index[id2] < 0 ) continue;
               alignlet la1 = aligns0[ aligns0_index[id1] ];
               alignlet la2 = aligns0[ aligns0_index[id2] ];
               int t1 = la1.TargetId( ), t2 = la2.TargetId( );
               if ( to_scaffold[t1] < 0 || to_scaffold[t2] < 0 ) continue;
               int sid = to_scaffold[t1];
               if ( to_scaffold[t2] != sid || !la1.Fw1( ) || la2.Fw1( ) ) continue;
               const superb& s = scaffolds[sid];
               int p1 = to_scaffold_pos[t1], p2 = to_scaffold_pos[t2];
               if ( !( p1 <= p2 ) ) continue;
               int start = la1.pos2( ), stop = la2.Pos2( );
               if ( p1 == p2 && !( start <= stop ) ) continue;
               int psep = pairs.sep( pairs.getPairID(id1) );
               int pdev = pairs.sd( pairs.getPairID(id1) );
               int sep, dev = pdev;
               if ( p1 == p2 ) sep = la2.pos2( ) - la1.Pos2( );
               else
               {    sep = s.SubSuperLength( p1, p2 ) - s.Len(p1) - s.Len(p2)
                         + la2.pos2( ) + tig_sizes[t1] - la1.Pos2( );
                    dev += s.SubSuperLengthDev( p1, p2 );    }
               if ( double( Abs( sep - psep ) ) > dev_mult * double(dev) ) continue;
               start += trim_back;
               stop -= trim_back;
               if ( p1 == p2 ) 
               {    if ( start <= stop ) covi[bi][t1].push( start, stop );    }
               else
               {    if ( start <= tig_sizes[t1] ) 
                         covi[bi][t1].push( start, tig_sizes[t1] );
                    for ( int j = p1 + 1; j < p2; j++ )
                         covi[bi][ s.Tig(j) ].push( 0, tig_sizes[ s.Tig(j) ] );
                    if ( stop >= 0 ) covi[bi][t2].push( 0, stop );    }    }    }
     cout << Date( ) << ": combining coverage" << endl;
     vec< vec<ho_interval> > cov(ntigs);
     #pragma omp parallel for
     for ( int t = 0; t < ntigs; t++ )
     {    for ( size_t bi = 0; bi < nbatches; bi++ )
               cov[t].append( covi[bi][t] );    }

     // Find uncovered parts of each edge.

     cout << Date( ) << ": finding uncovered parts" << endl;
     for ( int i = 0; i < trace_ids.isize( ); i++ )
     {    int t = trace_ids[i];
          Sort( cov[t] );
          log << "\ncoverage of " << t << endl;
          for ( int j = 0; j < cov[t].isize( ); j++ )
               log << cov[t][j] << "\n";    }
     vec< vec<ho_interval> > un(ntigs);
     #pragma omp parallel for
     for ( int t = 0; t < ntigs; t++ )
          Uncovered( tig_sizes[t], cov[t], un[t] );

     // Find suspicious gaps between contigs.  Not really handling tiny contigs
     // correctly.  This is compensating for a defect in scaffolding.

     cout << Date( ) << ": finding suspicious gaps" << endl;
     vec< vec<int> > suspicious_gaps( scaffolds.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < scaffolds.isize( ); i++ )
     {    const superb& s = scaffolds[i];
          for ( int j = 0; j < s.Ntigs( ); j++ )
          {    int t = s.Tig(j);
               if ( tig_sizes[t] < min_contig ) continue;
               if ( j < s.Ntigs( ) - 1 )
               {    if ( tig_sizes[ s.Tig(j+1) ] < min_contig ) continue;
                    vec<int> to_end;
                    for ( int l = 0; l < cov[t].isize( ); l++ )
                    {    if ( cov[t][l].Stop( ) == tig_sizes[t] )
                              to_end.push_back( cov[t][l].Start( ) );    }
                    Sort(to_end);
                    if ( to_end.empty( ) || to_end.back( ) - to_end.front( ) 
                         < min_spread )
                    {    log << "suspicious gap: " << t << " --> "
                              << s.Tig(j+1) << endl;    
                         suspicious_gaps[i].push_back(j);    }    }
               if ( j > 0 )
               {    if ( tig_sizes[ s.Tig(j-1) ] < min_contig ) continue;
                    vec<int> to_end;
                    for ( int l = 0; l < cov[t].isize( ); l++ )
                    {    if ( cov[t][l].Start( ) == 0 )
                              to_end.push_back( cov[t][l].Stop( ) );    }
                    Sort(to_end);
                    if ( to_end.empty( ) || to_end.back( ) - to_end.front( ) 
                         < min_spread )
                    {    log << "suspicious gap: " << s.Tig(j-1) << " --> "
                              << t << endl;    
                         suspicious_gaps[i].push_back(j-1);    }    }    }    }
     #pragma omp parallel for
     for ( int s = 0; s < scaffolds.isize( ); s++ )
          UniqueSort( suspicious_gaps[s] );

     // Find links that suggest misassambly.

     cout << Date( ) << ": finding evil links" << endl;
     vec< vec< vec<blink_fw> > > blink_fwsi(nbatches);
     vec< vec< vec<blink_rc> > > blink_rcsi(nbatches);
     for ( size_t j = 0; j < nbatches; j++ )
     {    blink_fwsi[j].resize(ntigs), blink_rcsi[j].resize(ntigs);    }
     #pragma omp parallel for
     for ( size_t bi = 0; bi < nbatches; bi++ )
     {    size_t start = ( bi * nreads ) / nbatches;
          size_t stop = ( (bi+1) * nreads ) / nbatches;
          for ( size_t id1 = start; id1 < stop; id1++ )
          {    if ( !pairs.isPaired(id1) ) continue;
               int id2 = pairs.getPartnerID(id1);
               if ( aligns0_index[id1] < 0 || aligns0_index[id2] < 0 ) continue;
               alignlet la1 = aligns0[ aligns0_index[id1] ];
               alignlet la2 = aligns0[ aligns0_index[id2] ];
               int t1 = la1.TargetId( ), t2 = la2.TargetId( );
	       ForceAssertLt( t1, ntigs );
	       ForceAssertLt( t2, ntigs );
               if ( to_scaffold[t1] < 0 || to_scaffold[t2] < 0 ) continue;
               int sid1 = to_scaffold[t1], sid2 = to_scaffold[t2];
               if ( sid1 == sid2 ) continue; // not really right
               const superb &s1 = scaffolds[sid1], &s2 = scaffolds[sid2];
               int p1 = to_scaffold_pos[t1], p2 = to_scaffold_pos[t2];
               int read_length2 = la2.Pos2( ) - la2.pos2( );
               int psep = pairs.sep( pairs.getPairID(id1) );
               int pdev = pairs.sd( pairs.getPairID(id1) );
     
               if ( la1.Fw1( ) )
               {    
                    // Ignore links that go off the end of the scaffold.
     
                    if ( la1.Pos2() + psep + int(round(dev_mult*pdev)) + read_length2
                         > s1.SubSuperLength( p1, s1.Ntigs( ) - 1 ) )
                    {    continue;    }
     
                    blink_fw b;
                    b.start1 = la1.pos2( );
                    b.stop1_low 
                         = b.start1 + psep - blink_dev_mult * pdev + read_length2;
                    b.stop1_high 
                         = b.start1 + psep + blink_dev_mult * pdev + read_length2;
                    b.t2 = t2;
                    b.stop2 = la2.Pos2( );
                    blink_fwsi[bi][t1].push_back(b);    }

               else
               {    
                    // Ignore links that go off the end of the scaffold.
     
                    if ( la1.pos2() - psep - int(round(dev_mult*pdev)) - read_length2
                         < 0 ) 
                    {    continue;    }
     
                    blink_rc b;
                    b.stop1 = la1.Pos2( );
                    b.start1_low 
                         = b.stop1 - psep - blink_dev_mult * pdev - read_length2;
                    b.start1_high 
                         = b.stop1 - psep + blink_dev_mult*pdev - read_length2;
                    b.t2 = t2;
                    b.start2 = la2.pos2( );
                    blink_rcsi[bi][t1].push_back(b);    }    }    }
     cout << Date( ) << ": combining evil links" << endl;
     vec< vec<blink_fw> > blink_fws(ntigs);
     vec< vec<blink_rc> > blink_rcs(ntigs);
     #pragma omp parallel for
     for ( int t = 0; t < ntigs; t++ )
     {    for ( size_t bi = 0; bi < nbatches; bi++ )
          {    blink_fws[t].append( blink_fwsi[bi][t] );
               blink_rcs[t].append( blink_rcsi[bi][t] );    }    }

     // Display uncovered parts of scaffolds and misassembly evidence associated
     // with them.  Find places where the scaffolds should be broken.

     cout << Date( ) << ": displaying uncovered parts" << endl;
     log << "\nuncovered parts of edges:\n";
     vec< vec<ho_interval> > splits(ntigs);
     for ( int i = 0; i < scaffolds.isize( ); i++ )
     {    const superb& s = scaffolds[i];
          if ( s.FullLength( ) < min_scaffold ) continue;
          log << "\nscaffold " << i << "\n";
          for ( int l = 0; l < s.Ntigs( ); l++)
          {    int t = s.Tig(l);
               log << t << "[l=" << tig_sizes[t] << "]:";
               for ( int j = 0; j < un[t].isize( ); j++ )
               {    const ho_interval& u = un[t][j];
                    if ( u.Stop( ) + s.SubSuperLength( 0, l ) - s.Len(l) 
                         < min_dist_from_end )
                    {    continue;    }
                    if ( -u.Start( ) + s.SubSuperLength( l, s.Ntigs( ) - 1 ) 
                         < min_dist_from_end )
                    {    continue;    }
                    log << " " << u;    

                    // Search for supporting links.
     
                    log << " {";
                    vec<blink_fw> fws;
                    for ( int l = 0; l < blink_fws[t].isize( ); l++ )
                    {    const blink_fw& b = blink_fws[t][l];
                         if ( !( b.start1 < u.Start( ) ) ) continue;
                         if ( !( b.stop1_high > u.Stop( ) ) ) continue;
                         fws.push_back(b);    }
                    Sort(fws);
                    Bool support = False;
                    for ( int l = 0; l < fws.isize( ); l++ )
                    {    int m;
                         for ( m = l + 1; m < fws.isize( ); m++ )
                              if ( fws[m].t2 != fws[l].t2 ) break;
                         int count = 0;
                         for ( int v = l; v < m; v++ )
                         {    if ( v == l || fws[v].stop2 != fws[v-1].stop2 )
                                   ++count;    }
                         if ( count < min_reach_away ) continue;
                         support = True;
                         log << " " << fws[l].t2 << "<";
                         for ( int v = l; v < m; v++ )
                         {    log << "(" << fws[v].start1 << ","
                                   << fws[v].stop2 << ")";    }
                         log << ">";
                         l = m - 1;    }
                    log << " ;";
                    vec<blink_rc> rcs;
                    for ( int l = 0; l < blink_rcs[t].isize( ); l++ )
                    {    const blink_rc& b = blink_rcs[t][l];
                         if ( !( b.stop1 > u.Stop( ) ) ) continue;
                         if ( !( b.start1_low < u.Start( ) ) ) continue;
                         rcs.push_back(b);    }
                    Sort(rcs);
                    for ( int l = 0; l < rcs.isize( ); l++ )
                    {    int m;
                         for ( m = l + 1; m < rcs.isize( ); m++ )
                              if ( rcs[m].t2 != rcs[l].t2 ) break;
                         int count = 0;
                         for ( int v = l; v < m; v++ )
                         {    if ( v == l || rcs[v].start2 != rcs[v-1].start2 )
                                   ++count;    }
                         if ( count < min_reach_away ) continue;
                         support = True;
                         log << " " << rcs[l].t2 << "<";
                         for ( int v = l; v < m; v++ )
                         {    log << "(" << rcs[v].stop1 << ","
                                   << rcs[v].start2 << ")";    }
                         log << ">";
                         l = m - 1;    }
                    log << " }";    
                    if (support) splits[t].push_back(u);    }
               log << "\n";    }    }

     // Process splits, yielding new scaffolds.
     
     int n_new_contigs = 0;
     
     cout << Date( ) << ": generating contigs and assembly fasta files" << endl;
     int scount = 0;
     for ( int si = 0; si < scaffolds.isize( ); si++ )
     {    const superb& S = scaffolds[si];
          if ( S.Ntigs() == 0 ) continue;
          vec<ho_interval> sp;
          int pos = 0;
          for ( int j = 0; j < S.Ntigs( ); j++ )
          {    int t = S.Tig(j);
               for ( int i = 0; i < splits[t].isize( ); i++ )
                    sp.push_back( splits[t][i] + pos );
               
               // Here we would cut at suspicious gaps, but we've turned this
               // code off.

               /*
               if ( Member( suspicious_gaps[si], j ) )
                    sp.push( pos, pos + Max( 1, S.Gap(j) ) );
               */

               if ( j < S.Ntigs( ) - 1 ) pos += S.Len(j) + Max( 1, S.Gap(j) );    }
          vec<char> s;
          // vector to indicate the beginning/end of the old contigs
          // on this scaffold
	  vec<int> contig_extents; 
          for ( int j = 0; j < S.Ntigs( ); j++ ) 
          {    fastavector b = contigs[ S.Tig(j) ];
	       for ( unsigned int l = 0; l < b.size( ); l++ ) 
               {    s.push_back( b[l] );
	            int extent_id = l == 0 ? (j+1) : l == b.size() - 1 ? -(j+1) : 0;
	            contig_extents.push_back( extent_id );    }
	       if ( j < S.Ntigs( ) - 1 ) 
               {    int gap = Max( 1, S.Gap(j) );
	            s.push_back_copies( 'n', gap );
	            contig_extents.push_back_copies( 0, gap );    }    }
          vec<ho_interval> un;
	  if ( sp.size() ) ForceAssertGt( s.isize(), sp.back().Stop() );
          Uncovered( s.size( ), sp, un );
	  
          // Make this scaffold, creating a superb object and printing the
          // bases to the fasta file.

          for ( int i = 0; i < un.isize( ); i++ )
          {
               // Remove gaps from the end of the scaffold.

               while( un[i].Start( ) < un[i].Stop( ) && s[ un[i].Start( ) ] == 'n' )
	            un[i].AddToStart(1);
               while( un[i].Start( ) < un[i].Stop( ) && s[ un[i].Stop( )-1 ] == 'n' )
	            un[i].AddToStop(-1);
               if ( un[i].Length( ) < min_scaffold ) continue;
               ForceAssertGe( (int) s.size(), un[i].Stop() );
               ForceAssertNe( s[ un[i].Start() ], 'n' );
               ForceAssertNe( s[ un[i].Stop()-1], 'n' );
       
               // Write to fasta (assembly).

               sout << ">scaffold_" << scount++ << "\n";
               int printed = 1;
               for ( int j = un[i].Start( ); j < un[i].Stop( ); j++, printed++ ) {
		 sout << s[j];
		 if ( printed % 80 == 0 ) sout << "\n";
	       }
               if ( printed == 1 || printed % 80 != 1 ) sout << endl;
	       
               // Break the scaffold on gaps ('n').

	       vec<int> cbegs;
	       vec<int> cends;
	       vec<int> gaps;    // has size cbegs.size( ) - 1
	       
	       int cbeg = un[i].Start( );
	       int pos = un[i].Start( );
	       while ( pos < un[i].Stop( ) ) {
		 while ( pos < un[i].Stop( ) && s[pos] != 'n' ) pos++;
		 cbegs.push_back( cbeg );
		 cends.push_back( pos );

		 int gap_begin = pos;
		 while ( pos < un[i].Stop( ) && s[pos] == 'n' ) pos++;
		 if ( pos != gap_begin ) gaps.push_back( pos - gap_begin );
		 cbeg = pos;
	       }
	       ForceAssertEq( cbegs.size( ), cends.size( ) );
	       ForceAssertEq( cbegs.size( ) - 1, gaps.size( ) );

               // Create the superb, and write to fasta (contigs).
	       
               superb new_scaffold;
               for ( int j = 0; j < cbegs.isize( ); j++ ) {
		 int clen = cends[j] - cbegs[j];;
		 int cstart = cbegs[j];
		 int cstop = cends[j]-1;
		 int gap = j > 0 ? gaps[j-1] : -1;
		 ForceAssertEq( clen, 1 + cstop - cstart );

		 cgout << ">contig_" << n_new_contigs << "\n";
		 int printed = 1;
		 for(int kk=cstart; kk<=cstop; kk++, printed++) {
		   cgout << s[kk];
		   if ( printed % 80 == 0 ) cgout << "\n";
		 }
		 if ( printed == 1 || printed % 80 != 1 ) cgout << "\n";
		 
		 if ( j == 0 ) new_scaffold.PlaceFirstTig( n_new_contigs++, clen );
		 else new_scaffold.AppendTig( n_new_contigs++, clen, gap, 0 );
	       }
               new_scaffolds.push_back( new_scaffold );
	  }
     }

}

