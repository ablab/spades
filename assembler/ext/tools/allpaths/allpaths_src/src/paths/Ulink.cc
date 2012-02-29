///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "paths/Ulink.h"

void DeleteBi( digraphE<ulink>& G, const int bi_mult, const Bool log )
{    vec<int> to_delete_bi;
     for ( int a = 0; a < G.N( ); a++ )
     {    for ( int j = 0; j < G.From(a).isize( ); j++ )
          {    int b = G.From(a)[j];
               const ulink& e1 = G.EdgeObjectByIndexFrom( a, j );
               for ( int l = 0; l < G.From(b).isize( ); l++ )
               {    if ( G.From(b)[l] != a ) continue;
                    const ulink& e2 = G.EdgeObjectByIndexFrom( b, l );
                    if ( e2.nlinks >= bi_mult * e1.nlinks )
                    {    to_delete_bi.push_back( 
                              G.EdgeObjectIndexByIndexFrom( a, j ) );
                         if (log)
                         {    cout << "deleting edge from " << a << " to " << b
                                   << " because the other direction is much "
                                   << "stronger\n";    }
                         goto to_delete_bi_next;    }    }
               to_delete_bi_next: continue;    }    }
     G.DeleteEdges(to_delete_bi);    }

void DeleteBiAnother( digraphE<ulink>& G, const vecbasevector& unibases,
     const Bool log )
{    vec<int> to_delete_bi;
     ulink nolink;
     for ( int a = 0; a < G.N( ); a++ )
     {    for ( int j = 0; j < G.From(a).isize( ); j++ )
          {    int b = G.From(a)[j];
               const ulink* ab = &G.EdgeObjectByIndexFrom( a, j );
               for ( int q = 0; q < G.From(b).isize( ); q++ )
               {    if ( G.From(b)[q] != a ) continue;
                    const ulink* ba = &G.EdgeObjectByIndexFrom( b, q );
                    int vote1 = 0, vote2 = 0;
                    vec<int> cs;
                    for ( int r = 0; r < G.To(a).isize( ); r++ )
                         cs.push_back( G.To(a)[r] );
                    for ( int r = 0; r < G.From(a).isize( ); r++ )
                         cs.push_back( G.From(a)[r] );
                    UniqueSort(cs);
                    int infinity = 1000000000;
                    double best_ab = infinity, best_ba = infinity;
                    for ( int l = 0; l < cs.isize( ); l++ )
                    {    int c = cs[l];
                         const ulink *ac = 0, *ca = 0, *bc = 0, *cb = 0;
                         for ( int s = 0; s < G.To(a).isize( ); s++ )
                         {    if ( G.To(a)[s] == c )
                              {    ca = &G.EdgeObjectByIndexTo( a, s );
                                   break;    }    }
                         for ( int s = 0; s < G.From(a).isize( ); s++ )
                         {    if ( G.From(a)[s] == c )
                              {    ac = &G.EdgeObjectByIndexFrom( a, s );
                                   break;    }    }
                         for ( int s = 0; s < G.From(b).isize( ); s++ )
                         {    if ( G.From(b)[s] == c )
                              {    bc = &G.EdgeObjectByIndexFrom( b, s );
                                   break;    }    }
                         for ( int s = 0; s < G.To(b).isize( ); s++ )
                         {    if ( G.To(b)[s] == c )
                              {    cb = &G.EdgeObjectByIndexTo( b, s );
                                   break;    }    }
                         // 111: a --> b, a --> c, b --> c
                         if ( ac != 0 && bc != 0 )
                         {    double err = ac->Sep( ) - ab->Sep( ) - bc->Sep( )
                                   - unibases[b].isize( );
                              double dev = sqrt( double( ac->Var( ) + ab->Var( )
                                   + bc->Var( ) ) );
                              best_ab = Min( best_ab, Abs(err/dev) );    }
                         // 112: a --> b, a --> c, c --> b
                         if ( ac != 0 && cb != 0 )
                         {    double err = ab->Sep( ) - ac->Sep( ) - cb->Sep( )
                                   - unibases[c].isize( );
                              double dev = sqrt( double( ab->Var( ) + ac->Var( )
                                   + cb->Var( ) ) );
                              best_ab = Min( best_ab, Abs(err/dev) );    }
                         // 122: a --> b, c --> a, c --> b
                         if ( ca != 0 && cb != 0 )
                         {    double err = cb->Sep( ) - ca->Sep( ) - ab->Sep( )
                                   - unibases[a].isize( );
                              double dev = sqrt( double( cb->Var( ) + ca->Var( )
                                   + ab->Var( ) ) );
                              best_ab = Min( best_ab, Abs(err/dev) );    }
                         // 211: b --> a, a --> c, b --> c
                         if ( ac != 0 && bc != 0 )
                         {    double err = bc->Sep( ) - ba->Sep( ) - ac->Sep( )
                                   - unibases[a].isize( );
                              double dev = sqrt( double( bc->Var( ) + ba->Var( )
                                   + ac->Var( ) ) );
                              best_ba = Min( best_ba, Abs(err/dev) );    }
                         // 221: b --> a, c --> a, b --> c
                         if ( ca != 0 && bc != 0 )
                         {    double err = ba->Sep( ) - bc->Sep( ) - ca->Sep( )
                                   - unibases[c].isize( );
                              double dev = sqrt( double( ba->Var( ) + bc->Var( )
                                   + ca->Var( ) ) );
                              best_ba = Min( best_ba, Abs(err/dev) );    }
                         // 222: b --> a, c --> a, c --> b
                         if ( ca != 0 && cb != 0 )
                         {    double err = ca->Sep( ) - cb->Sep( ) - ba->Sep( )
                                   - unibases[b].isize( );
                              double dev = sqrt( double( ca->Var( ) + cb->Var( )
                                   + ba->Var( ) ) );
                              best_ba = Min( best_ba, Abs(err/dev) );    }
                         if ( best_ab + 3.0 <= best_ba ) vote1++;
                         if ( best_ba + 3.0 <= best_ab ) vote2++;
                         if ( log && ( best_ab < infinity || best_ba < infinity ) )
                         {    cout << "DeleteBiAnother evidence: ";
                              PRINT5( a, b, c, best_ab, best_ba );    }    }
                    if ( vote1 > 0 && vote2 == 0 )
                    {    to_delete_bi.push_back( 
                              G.EdgeObjectIndexByIndexFrom( b, q ) );    
                         if (log)
                         {    cout << "DeleteBiAnother: deleting edge from "
                                   << b << " to " << a << "\n";    }    }
                    if ( vote2 > 0 && vote1 == 0 )
                    {    to_delete_bi.push_back( 
                              G.EdgeObjectIndexByIndexFrom( a, j ) );
                         if (log)
                         {    cout << "DeleteBiAnother: deleting edge from " << a 
                                   << " to " << b << "\n";    }    }    }    }    }
     G.DeleteEdges(to_delete_bi);    }

void DeleteNlinks( digraphE<ulink>& G, const int min_links, const Bool log )
{    vec<int> to_delete_nlinks;
     for ( int b = 0; b < G.N( ); b++ )
     {    for ( int j = 0; j < G.From(b).isize( ); j++ )
          {    int c = G.From(b)[j];
               const ulink& e2 = G.EdgeObjectByIndexFrom( b, j );
               if ( e2.nlinks >= min_links ) continue;
               Bool accept = False;
               for ( int l = 0; l < G.To(b).isize( ); l++ )
               {    int a = G.To(b)[l];
                    const ulink& e1 = G.EdgeObjectByIndexTo( b, l );
                    if ( e1.nlinks < min_links ) continue;
                    for ( int m = 0; m < G.From(a).isize( ); m++ )
                    {    if ( G.From(a)[m] != c ) continue;
                         const ulink& e3 = G.EdgeObjectByIndexFrom( a, m );
                         if ( e3.nlinks < min_links ) continue;
                         accept = True;
                         goto nlinks_tail;    }    }
               nlinks_tail: if (accept) continue;
               to_delete_nlinks.push_back( 
                    G.EdgeObjectIndexByIndexFrom( b, j ) );    }    }
     G.DeleteEdges(to_delete_nlinks);    }

void DeleteTrans( digraphE<ulink>& G, const vecbasevector& unibases, 
     const int trans_depth, const int max_sep, const Bool log )
{    vec<int> to_delete_trans;
     for ( int v = 0; v < G.N( ); v++ )
     {    for ( int j = 0; j < G.From(v).isize( ); j++ )
          {    int t = G.From(v)[j];
               const ulink& e_trans = G.EdgeObjectByIndexFrom( v, j );
               if ( e_trans.Sep( ) - e_trans.Dev( ) < max_sep ) continue;
               vec< vec<int> > paths;
               for ( int k = 0; k < G.From(v).isize( ); k++ )
               {    int w = G.From(v)[k];
                    if ( w != t ) 
                    {    vec<int> p(2);
                         p[0] = v, p[1] = w;
                         paths.push_back(p);    }    }
               int k = 0;
               for ( int tpass = 1; tpass < trans_depth; tpass++ )
               {    int npaths = paths.size( );
                    for ( ; k < npaths; k++ )
                    {    vec<int> p = paths[k];
                         if ( paths[k].back( ) != t )
                         {    int x = paths[k].back( );
                              for ( int m = 0; m < G.From(x).isize( ); m++ )
                              {    vec<int> p = paths[k];
                                   p.push_back( G.From(x)[m] );
                                   paths.push_back(p);    }    }    }    }
               Bool trans = False;
               for ( int k = 0; k < paths.isize( ); k++ )
                    if ( paths[k].back( ) == t ) trans = True;
               if ( !trans ) continue;
               to_delete_trans.push_back( G.EdgeObjectIndexByIndexFrom( v, j ) );
               vec< pair<int,int> > acceptors;
               for ( int k = 0; k < paths.isize( ); k++ )
               {    if ( paths[k].back( ) == t )
                    {    for ( int l = 1; l < paths[k].isize( ); l++ )
                              acceptors.push( paths[k][l-1], paths[k][l] );    }    }
               UniqueSort(acceptors);
               if (log)
               {    cout << "deleting link from " << v << " to " << t
                         << " based on transitivity, using";
                    for ( int k = 0; k < acceptors.isize( ); k++ )
                    {    int x = acceptors[k].first, y = acceptors[k].second;
                         cout << " (" << x << "," << y << ")";    }
                    cout << "\n";    }
               for ( int k = 0; k < acceptors.isize( ); k++ )
               {    int x = acceptors[k].first, y = acceptors[k].second;
                    vec<int> es = G.EdgesBetween( x, y );
                    ForceAssertEq( es.isize( ), 1 );
                    ulink& e = G.EdgeObjectMutable(es[0]);
                    if ( x == v ) e.start1 = Min( e.start1, e_trans.start1 );
                    else e.start1 = 0;
                    if ( y == t ) e.stop2 = Max( e.stop2, e_trans.stop2 );
                    else e.stop2 = unibases[y].size( );    }    }    }
     G.DeleteEdges(to_delete_trans);    }

void DeleteWeak( digraphE<ulink>& G, const vecbasevector& unibases,
     const int min_glue, const Bool log )
{    vec<int> to_delete_weak;
     for ( int u1 = 0; u1 < G.N( ); u1++ )
     {    for ( int j = 0; j < G.From(u1).isize( ); j++ )
          {    int u2 = G.From(u1)[j];
               const ulink& e = G.EdgeObjectByIndexFrom( u1, j );
               Bool good_left = ( e.start1 == 0 
                    || unibases[u1].isize( ) - e.start1 >= min_glue );
               Bool good_right = ( e.stop2 == unibases[u2].isize( )
                    || e.stop2 >= min_glue );
               if ( !good_left || !good_right )
               {    if (log)
                    {    cout << "deleting weak link from " << u1
                              << " to " << u2 << "\n";    }
                    to_delete_weak.push_back( 
                         G.EdgeObjectIndexByIndexFrom( u1, j ) );    }    }    }
     G.DeleteEdges(to_delete_weak);    }
