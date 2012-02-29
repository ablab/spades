///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/BigMapTools.h"
#include "paths/Uniseq.h"

void BigMapDot( const String& DOT, const snark& S, const int K2, 
     const Bool VALIDATE, const vecbasevector& genome, 
     const vec< vec<placementy> >& Ulocs2, const vec<String>& legends_to_show,
     const Bool circo )
{    
     // Set up.

     vec< vec<String> > legends;
     vec<String> colors;

     // Generate legend.

     vec<String> vertex_labels, legend1;
     legend1.push_back( "<FONT POINT-SIZE=\"12\">Legend</FONT>" );
     int count = 0;
     for ( int x = 0; x < S.VertN( ); x++ )
     {    if ( S.Vert(x).N( ) == 1 )
               vertex_labels.push_back( ToString( S.Vert(x).U(0) ) );
          else
          {    String s = BaseAlpha(count) + " = ";
               int page_width = 120;
               for ( int j = 0; j < S.Vert(x).U( ).isize( ); j++ )
               {    if ( j > 0 ) 
                    {    int over = S.Vert(x).Over(j-1);
                         if ( over == K2 - 1 ) s += "+";
                         else s += "--(" + ToString(over) + ")--";     }
                    s += ToString( S.Vert(x).U(j) );
                    if ( j < S.Vert(x).N( ) - 1 && s.isize( ) >= page_width )
                    {    legend1.push_back(s);
                         s = "    ";    }    }
               vertex_labels.push_back( BaseAlpha(count++) );
               legend1.push_back(s);    }    }
     vec< vec<String> > edge_labels( S.VertN( ) );
     for ( int v = 0; v < S.VertN( ); v++ )
     {    for ( int j = 0; j < S.G( ).From(v).isize( ); j++ )
          {   int n = S.G( ).EdgeObjectByIndexFrom( v, j ).ClosureCount( );
               edge_labels[v].push_back( ToString(n) );    }    }

     Ofstream( dout, DOT );
     if ( Member( legends_to_show, String("Legend") ) )
     {    legends.push_back(legend1);
          colors.push_back( "aquamarine" );    }

     // Look for problems.

     int prob_count = 0, severe_prob_count = 0;
     if (VALIDATE)
     {    
          // Set up.

          vec<String> legend2;
          legend2.push_back( "<FONT POINT-SIZE=\"12\">Discrepancies</FONT>" );
          String circ = "&#8226; ";
          cout << "\n============================================================="
               << "=======================\n";
          cout << "\nPROBLEMS\n";

          // Find the unipaths that have to be used, either because they are in a
          // vertex, or because they are in the intersection of an edge object.

          vec<Bool> required( S.Unibases( ).size( ), False );
          for ( int x = 0; x < S.VertN( ); x++ )
          {    for ( int j = 0; j < S.Vert(x).N( ); j++ )
                    required[ S.Vert(x).U(j) ] = True;    }
          for ( int x = 0; x < S.EdgeN( ); x++ )
          {    if ( S.Edge(x).Closed( ) )
               {    vec<int> intersection = Common( S.Edge(x).Closures( ) );
                    for ( int j = 0; j < intersection.isize( ); j++ )
                         required[ intersection[j] ] = True;    }    }

          // Note problems for required unipaths that have errors (ignoring multiple 
          // placements).  Note that possibility of double counting if an error
          // occurs in the overlap between two unibases.

          const int infinity = 1000000000;
          vec<int> known_errors( S.Unibases( ).size( ), 0 );
          for ( int u = 0; u < required.isize( ); u++ )
          {    if ( !required[u] ) continue;
               if ( Ulocs2[u].size( ) >= 2 ) continue;
               int min_score = infinity;
               if ( Ulocs2[u].solo( ) ) min_score = Ulocs2[u][0].Errs( );
               if ( Ulocs2[u].empty( ) )
               {    cout << "\nPROBLEM " << ++prob_count << "\n";
                    cout << "from: " << u << "\n";
                    PRINT(min_score);
                    legend2.push_back( ToString(prob_count) + ". " 
                         + "{" + ToString(u) + "}, unmapped, errs = &infin;" );
                    known_errors[u] = min_score;    }
               else if ( min_score > 0 )
               {    cout << "\nPROBLEM " << ++prob_count << "\n";
                    cout << "from: " << u << "\n";
                    PRINT(min_score);
                    int mismatches = Ulocs2[u][0].mismatches;
                    int indels = Ulocs2[u][0].indels;
                    int unmapped = Ulocs2[u][0].unmapped;
                    String errs;
                    Bool printed = False;
                    if ( mismatches > 0 )
                    {    errs = ToString(mismatches) + " mismatches";
                         printed = True;    }
                    if ( indels > 0 )
                    {    if (printed) errs += " + ";
                         errs += ToString(indels) + " indels";
                         printed = True;    }
                    if ( unmapped > 0 )
                    {    if (printed) errs += " + ";
                         errs += ToString(unmapped) + " unmapped";    }
                    legend2.push_back( ToString(prob_count) + ". " 
                         + "{" + ToString(u) + "}, errs = " + errs );
                    known_errors[u] = min_score;    }    }

          // Find certain lines in the graph that have to map to somewhere on the
          // genome.  We start at a given edge, then go forward as far as possible,
          // then backward as far as possible.  A line is a list of consecutive
          // vertices.

          vec< vec<int> > lines;
          for ( int x = 0; x < S.VertN( ); x++ )
          for ( int j = 0; j < S.G( ).From(x).isize( ); j++ )
          {    vec<int> line;
               int e = S.G( ).EdgeObjectIndexByIndexFrom( x, j );
               int y = S.G( ).From(x)[j];
               line.push_back( x, y );
               Bool cycled = False;
               while( S.G( ).From(y).solo( ) )
               {    int e = S.G( ).EdgeObjectIndexByIndexFrom( y, 0 );
                    y = S.G( ).From(y)[0];
                    cycled = Member( line, y );
                    line.push_back(y);
                    if (cycled) break;    }
               if ( !cycled )
               {    y = x;
                    while( S.G( ).To(y).solo( ) )
                    {    int e = S.G( ).EdgeObjectIndexByIndexTo( y, 0 );
                         y = S.G( ).To(y)[0];
                         cycled = Member( line, y );
                         line.push_front(y);
                         if (cycled) break;    }    }
               lines.push_back(line);    }
          UniqueSort(lines);

          // Remove lines that are subsets of other lines.

          vec<Bool> to_delete( lines.size( ), False );
          for ( int i1 = 0; i1 < lines.isize( ); i1++ )
          for ( int i2 = 0; i2 < lines.isize( ); i2++ )
          {    if ( i2 == i1 ) continue;
               if ( lines[i1].Contains( lines[i2] ) ) to_delete[i2] = True;    }
          EraseIf( lines, to_delete );

          // Now translate to uniseqs.

          vec< vec< vec<uniseq> > > ulines( lines.size( ) );
          for ( int i = 0; i < lines.isize( ); i++ )
          {    vec< vec<uniseq> >& ul = ulines[i];
               for ( int j = 0; j < lines[i].isize( ); j++ )
               {    vec<uniseq> vert;
                    vert.push_back( S.Vert( lines[i][j] ) );
                    ulines[i].push_back(vert);
                    if ( j < lines[i].isize( ) - 1 )
                    {    int x1 = lines[i][j], x2 = lines[i][j+1];
                         for ( int k = 0; k < S.G( ).From(x1).isize( ); k++ )
                         {    if ( S.G( ).From(x1)[k] == x2 )
                              {    const gapster& g
                                        = S.G( ).EdgeObjectByIndexFrom( x1, k );
                                   ulines[i].push_back( g.Closures( ) );
                                   break;    }    }    }    }    }

          // Now split ulines wherever there is a uniquely placed unipath on a
          // vertex.

          vec< vec< vec<uniseq> > > ulines2;
          for ( int i = 0; i < ulines.isize( ); i++ )
          {    vec< vec<uniseq> >& ul = ulines[i];
               restart:
               for ( int j = 0; j < ul.isize( ); j += 2 )
               {    const uniseq& us = ul[j][0];
                    for ( int k = 0; k < us.N( ); k++ )
                    {    int u = us.U(k);
                         if ( !Ulocs2[u].solo( ) ) continue;
                         if ( k == 0 && j == 0 ) continue;
                         if ( k == us.N( ) - 1 && j == ul.isize( ) - 1 ) continue;
                         vec< vec<uniseq> > head, tail;
                         for ( int l = 0; l < j; l++ )
                              head.push_back( ul[l] );
                         vec<int> uh = us.U( ), oh = us.Over( );
                         uh.resize(k+1), oh.resize(k);
                         uniseq uhead( uh, oh );
                         vec<int> ut, ot;
                         ut.SetToSubOf( us.U( ), k, us.N( ) - k );
                         ot.SetToSubOf( us.Over( ), k, us.N( ) - k - 1 );
                         uniseq utail( ut, ot );
                         vec<uniseq> vuhead, vutail;
                         vuhead.push_back(uhead), vutail.push_back(utail);
                         head.push_back(vuhead), tail.push_back(vutail);
                         for ( int l = j+1; l < ul.isize( ); l++ )
                              tail.push_back( ul[l] );
                         ulines2.push_back(head);
                         ul = tail;
                         goto restart;    }    }
               ulines2.push_back(ul);    }
          ulines = ulines2;
          UniqueSort(ulines);

          // Remove subset ulines.

          to_delete.resize_and_set( ulines.size( ), False );
          for ( int i1 = 0; i1 < ulines.isize( ); i1++ )
          for ( int i2 = 0; i2 < ulines.isize( ); i2++ )
          {    if ( i2 == i1 ) continue;
               if ( ulines[i1].Contains( ulines[i2] ) ) to_delete[i2] = True;    }
          EraseIf( ulines, to_delete );

          // Now expand out the ulines, so that each becomes a vec of
          // uniseq --gap--> uniseq --gap--> ....  This is exponential, so could
          // be a problem.  The new data structure is also a 
          // vec< vec< vec<uniseq> > >, but it means something different.

          vec< vec< vec<uniseq> > > zlines;
          for ( int i = 0; i < ulines.isize( ); i++ )
          {    const vec< vec<uniseq> >& ul = ulines[i];

               // Do a quick precheck on the size of the structure we're going to
               // build.  If it's too big, exit ungracefully.

               const int64_t max_children = 1000000;
               int64_t M = 1;
               for ( int j = 0; j < ul.isize( ); j++ )
               {    M *= ul[j].size( );
                    if ( M > max_children ) break;    }
               if ( M > max_children )
               {    cout << "\nBigMapDot: number of children exceeded, giving up."
                         << endl;
                    _exit(1);    }
               vec< vec<uniseq> > zl;
               zl.push_back( ul[0] );
               for ( int j = 1; j < ul.isize( ); j++ )
               {    vec< vec<uniseq> > zlnew;
                    if ( ul[j].empty( ) ) continue;
                    else if ( ul[j-1].nonempty( ) )
                    {    for ( int k = 0; k < zl.isize( ); k++ )
                         {    for ( int l = 0; l < ul[j].isize( ); l++ )
                              {    vec<uniseq> z = zl[k];
                                   z.back( ) = Cat( z.back( ), ul[j][l] );
                                   zlnew.push_back(z);    }    }    }
                    else
                    {    for ( int k = 0; k < zl.isize( ); k++ )
                         {    for ( int l = 0; l < ul[j].isize( ); l++ )
                              {    vec<uniseq> z = zl[k];
                                   z.push_back( ul[j][l] );
                                   zlnew.push_back(z);    }    }    }
                    zl = zlnew;    }
               zlines.push_back(zl);    }
          for ( int i = 0; i < zlines.isize( ); i++ )
          {    
               // OK, is it any good?

               const int max_placements = 100000;
               int min_score = infinity;
               for ( int j = 0; j < zlines[i].isize( ); j++ )
               {    const vec<uniseq>& z = zlines[i][j];
                    vec< vec<placementy> > placements;
                    int u0 = z[0].U(0);
                    placements.push_back( vec<placementy>( ) );
                    for ( int k = 0; k < z.isize( ); k++ )
                    {    for ( int l = 0; l < z[k].N( ); l++ )
                         {    vec< vec<placementy> > placements_new;
                              int u1 = z[k].U(l);
                              for ( int r = 0; r < Ulocs2[u1].isize( ); r++ )
                              {    const placementy& p2 = Ulocs2[u1][r];
                                   #pragma omp parallel for
                                   for ( int s = 0; s < placements.isize( ); s++ )
                                   {    vec<placementy> P = placements[s];
                                        if ( P.nonempty( ) )
                                        {    const placementy& p1 = P.back( );
                                             if ( p1.g != p2.g ) continue;
                                             if ( p1.fw != p2.fw ) continue;
                                             int sep1 = p1.fw ? p2.pos - p1.Pos 
                                                  : p1.pos - p2.Pos;
                                             int sep2;
                                             if ( p1.fw ) 
                                             {    sep2 = -genome[p1.g].isize( ) 
                                                       - p1.Pos + p2.pos;    }
                                             else 
                                             {    sep2 = -genome[p1.g].isize( ) 
                                                       - p2.Pos + p1.pos;    }
                                             int msep = Min( Abs(sep1), Abs(sep2) );
                                             // BAD, NOTE HARDCODED VALUES!
                                             if ( msep > 15000 ) continue;
                                             if ( l > 0 && msep > 1000 )
                                                  continue;    
                                                       }
                                        P.push_back(p2);
                                        #pragma omp critical
                                        placements_new.push_back(P);    }    }
                              placements = placements_new;    
                              if ( placements.isize( ) > max_placements )
                                   goto next_i;    }    }
                    vec<int> score( placements.size( ), 0 );
                    vec<Bool> invalid( placements.size( ), False );
                    ParallelSort(placements);
                    #pragma omp parallel for
                    for ( int m = 0; m < placements.isize( ); m++ )
                    {    int count = 0;
                         for ( int k = 0; k < z.isize( ); k++ )
                         {    for ( int l = 0; l < z[k].N( ); l++ )
                              {    const placementy& p2 = placements[m][count];
                                   score[m] += p2.Errs( );
                                   if ( count > 0 )
                                   {    const placementy& 
                                             p1 = placements[m][count-1];
                                        int sep1 = p1.fw ? p2.pos - p1.Pos 
                                             : p1.pos - p2.Pos;
                                        int sep2;
                                        if ( p1.fw ) 
                                        {    sep2 = -genome[p1.g].isize( ) 
                                                  - p1.Pos + p2.pos;    }
                                        else 
                                        {    sep2 = -genome[p1.g].isize( ) 
                                                  - p2.Pos + p1.pos;    }
                                        if ( l == 0 )
                                        {    // BAD BAD BAD, NOTE HARDCODED VALUES!
                                             if ( (sep1 < -K2+1 || sep1 > 15000)
                                                  && (sep2 < -K2+1 || sep2 > 15000) )
                                             {    invalid[m] = True;    }    }
                                        else
                                        {    int expected_sep = -z[k].Over(l-1);
                                             score[m] += Min( Abs(sep1-expected_sep),
                                                  Abs(sep2-expected_sep) );    }    }
                                   count++;    }    }    }
                    for ( int m = 0; m < placements.isize( ); m++ )
                         min_score = Min( min_score, score[m] );    }
               if ( min_score > 0 )
               {    int eknown = 0;
                    for ( int j = 0; j < ulines[i].isize( ); j++ )
                    {    vec< vec<int> > core;
                         for ( int k = 0; k < ulines[i][j].isize( ); k++ )
                         {    vec<int> x = ulines[i][j][k].U( );
                              if ( k % 2 == 1 ) 
                              {
                                   ForceAssertGe( x.isize( ), 2 ); // XXXXXXXXXXXXXX
                                   x = SubOf( x, 1, x.isize( ) - 2 );
                                        }
                              core.push_back(x);   }
                         if ( core.nonempty( ) )
                         {    vec<int> intersection;
                              Intersection( core, intersection );
                              for ( int k = 0; k < intersection.isize( ); k++ )
                              {    eknown += known_errors[ 
                                        intersection[k] ];    }    }    }
                    if ( eknown == min_score ) continue;
                    cout << "\nPROBLEM " << ++prob_count << "\n";
                    // PRINT( zlines[i].size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    if ( min_score == infinity ) severe_prob_count++;

                    /*
                    for ( int j = 0; j < zlines[i].isize( ); j++ ) // XXXXXXXXXXXXXX
                    {    cout << "[" << j+1 << "] "; // XXXXXXXXXXXXXXXXXXXXXXXXXXXX
                         for ( int k = 0; k < zlines[i][j].isize( ); k++ ) // XXXXXX
                         {    if ( k > 0 ) cout << " --> "; // XXXXXXXXXXXXXXXXXXXXX
                              zlines[i][j][k].Print( cout, K2 );    } // XXXXXXXXXXX
                         cout << endl;    } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    */

                    cout << "from: ";
                    Print( cout, ulines[i], K2 );
                    cout << endl;
                    PRINT(min_score);    

                    ostringstream uout;
                    Print( uout, ulines[i], K2 );
                    String s = uout.str( );
                    const int pw = 60;
                    if ( s.isize( ) > pw )
                    {    s.resize(pw);
                         s += "...";    }
                    s.GlobalReplaceBy( "-->", "&rarr;" );
                    String errs = ( min_score < infinity ? ToString(min_score)
                         : "&infin;" );
                    legend2.push_back( ToString(prob_count) + ". " 
                         + s + ", errs = " + errs );    }
               next_i: continue;    }

          // Finish problem report.

          if ( Member( legends_to_show, String("Discrepancies") ) )
          {    legends.push_back(legend2);
               colors.push_back( "cyan" );    }    }

     // Generate summary.

     vec<String> legend3;
     legend3.push_back( "<FONT POINT-SIZE=\"12\">Summary</FONT>" );
     double genome_Mb = double( S.EstimatedGenomeSize( ) ) / 1000000.0;
     legend3.push_back( ToString( genome_Mb, 1 ) + " Mb" );
     legend3.push_back( ToString( S.EdgeN( ) ) + " edges" );
     int open = 0, ambig = 0;
     for ( int i = 0; i < S.EdgeN( ); i++ )
     {    const gapster& g = S.Edge(i);
          if ( g.Open( ) ) open++;
          else ambig += g.ClosureCount( ) - 1;    }
     legend3.push_back( ToString(open) + " gaps" );
     legend3.push_back( ToString(ambig) + " ambiguities" );
     if (VALIDATE) legend3.push_back( ToString(prob_count) + " discrepancies" );
     if (VALIDATE) 
          legend3.push_back( ToString(severe_prob_count) + " severe discrepancies" );
     if ( Member( legends_to_show, String("Summary") ) )
     {    legends.push_back(legend3);
          colors.push_back( "lightsalmon" );    }

     // Print dot:

     S.G( ).DOT_vl( dout, vertex_labels, ( circo ? "circo" : "" ), legends, colors, 
          edge_labels );    }
