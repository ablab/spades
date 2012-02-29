// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#include "ContigEvent.h"
#include "Set.h"
#include "String.h"
#include "Vec.h"

vec<String> contig_event::class_names;
vec<void*> contig_event::virtual_table_pointers;

#define declare_contig_event(EVENT_NAME)                                          \
{    contig_event::class_names.push_back( #EVENT_NAME );                          \
     EVENT_NAME ev;                                                               \
     contig_event::virtual_table_pointers.push_back( ( (void**) &ev )[0] );    }

void InitializeContigEventStaticMembers( )
{    if ( contig_event::class_names.empty( ) )
     {    declare_contig_event( reverse_contig_event );
          declare_contig_event( rename_contig_event );    
          declare_contig_event( merge_contig_event );
          declare_contig_event( merge_via_chimera_contig_event );
          declare_contig_event( delete_contig_event );
          declare_contig_event( link_contig_event );
          declare_contig_event( reverse_link_contig_event );
          declare_contig_event( mark_mergecontigs_start_contig_event );
          declare_contig_event( mark_mergesupercontigs_start_contig_event );   }   }

int contig_event::ContigEventClassId( ) const
{    InitializeContigEventStaticMembers( );
     String class_name = ContigEventClassName( );
     for ( unsigned int i = 0; i < class_names.size( ); i++ )
          if ( class_names[i] == class_name ) return i;
     FatalErr( "Please add the line\n"
          << "          declare_contig_event( " << class_name << " );\n" 
          << "to the definition of InitializeStaticMembers in ContigEvent.cc." );
     return -1;    }

void contig_event_list::Dump( )
{    for ( unsigned int i = 0; i < event_locs_.size( ); i++ )
          (*this)[i].Print(cout);    }

void contig_event_list::TraceBack( int e )
{    set<int> sources;
     sources.insert(e);
     vec<int> prints;
     for ( int i = (int) event_locs_.size( ) - 1; i >= 0; i-- )
     {    const contig_event& ce = (*this)[i];
          if ( ce.Marker( ) ) prints.push_back(i);
          else if ( ce.ContigOperator( ) )
          {    for ( set<int>::iterator j = sources.begin( ); j != sources.end( );
                    j++ )
                    if ( ce.IsOutputContig( *j ) )
                    {    prints.push_back(i);
                         vec<int> v = ce.InputContigs( );
                         sources.erase( *j );
                         for ( unsigned int k = 0; k < v.size( ); k++ )
                              sources.insert( v[k] );    
                         break;    }    }    }
     for ( int i = (int) prints.size( ) - 1; i >= 0; i-- )
          (*this)[ prints[i] ].Print(cout);    }

void contig_event_list::TraceBack( int e1, int e2 )
{    set<int> sources1, sources2;
     sources1.insert(e1);
     sources2.insert(e2);
     vec<int> prints;
     for ( int i = (int) event_locs_.size( ) - 1; i >= 0; i-- )
     {    const contig_event& ce = (*this)[i];
          if ( ce.Marker( ) ) prints.push_back(i);
          else if ( ce.ContigOperator( ) )
          {    for ( set<int>::iterator j = sources1.begin( ); j != sources1.end( );
                    j++ )
                    if ( ce.IsOutputContig( *j ) )
                    {    prints.push_back(i);
                         vec<int> v = ce.InputContigs( );
                         sources1.erase( *j );
                         for ( unsigned int k = 0; k < v.size( ); k++ )
                              sources1.insert( v[k] );    
                         break;    }
               for ( set<int>::iterator j = sources2.begin( ); j != sources2.end( );
                    j++ )
                    if ( ce.IsOutputContig( *j ) )
                    {    prints.push_back(i);
                         vec<int> v = ce.InputContigs( );
                         sources2.erase( *j );
                         for ( unsigned int k = 0; k < v.size( ); k++ )
                              sources2.insert( v[k] );    
                         break;    }    }
          else if ( ce.ReverseLinkOperator( ) )
          {    vec<int> v = ce.InputContigs( );
               if ( Member( sources1, v[1] ) && Member( sources2, v[0] ) )
               {    prints.push_back(i);
                    swap( sources1, sources2 );    }    }
          else if ( ce.LinkCreationOperator( ) )
          {    vec<int> v = ce.InputContigs( );
               if ( Member( sources1, v[0] ) && Member( sources2, v[1] ) )
                    prints.push_back(i);    }    }
     for ( int i = (int) prints.size( ) - 1; i >= 0; i-- )
          (*this)[ prints[i] ].Print(cout);    }
