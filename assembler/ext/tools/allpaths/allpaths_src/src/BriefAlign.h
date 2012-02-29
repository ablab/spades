// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#include <iostream>

#include "Misc.h"

#ifndef BRIEF_ALIGN
#define BRIEF_ALIGN

// Class brief_align encodes the minimal information regarding an alignment between
// two sequences.

class brief_align {

     public:

     int id1, id2;
     Bool rc2;
     int offset;

     brief_align( ) { }

     brief_align( int id1_arg,
		  int id2_arg,
		  Bool rc2_arg,
		  int offset_arg )
       : id1(id1_arg),
	 id2(id2_arg),
	 rc2(rc2_arg),
	 offset(offset_arg)
     { }

     int Id1( ) const { return id1; }
     int Id2( ) const { return id2; }
     Bool Rc2( ) const { return rc2; }
     int Offset( ) const { return offset; }

     friend Bool operator<( const brief_align& a1, const brief_align& a2 )
     {    return a1.id1 < a2.id1 
               || ( a1.id1 == a2.id1 
               && ( a1.id2 < a2.id2
               || ( a1.id2 == a2.id2 
               && ( a1.rc2 < a2.rc2
               || ( a1.rc2 == a2.rc2 && a1.offset < a2.offset ) ) ) ) );    }

     friend istream& operator>>( istream& i, brief_align& a )
     {    i.read( (char*) &a.id1, sizeof(a.id1) );
          i.read( (char*) &a.id2, sizeof(a.id2) );
          i.read( (char*) &a.rc2, sizeof(a.rc2) );
          i.read( (char*) &a.offset, sizeof(a.offset) );
          return i;    }

     friend ostream& operator<<( ostream& o, const brief_align& a )
     {    o.write( (char*) &a.id1, sizeof(a.id1) );
          o.write( (char*) &a.id2, sizeof(a.id2) );
          o.write( (char*) &a.rc2, sizeof(a.rc2) );
          o.write( (char*) &a.offset, sizeof(a.offset) );
          return o;    }

};

// This functor is useful for searching algorithms like equal_range.

#include <functional>

struct order_brief_align_by_id1
  : public binary_function<brief_align,brief_align,bool>
{
  bool operator() ( const brief_align& a1, const brief_align& a2) const
     {    
	 return a1.id1 < a2.id1 
               || ( a1.id1 == a2.id1 
	       && ( a1.id2 < a2.id2 
	       || ( a1.id2 == a2.id2 
	       && ( a1.rc2 < a2.rc2 ))));
     }
};

#endif
