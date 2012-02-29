///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Class map_int_int represents a function from the integers to the integers,
// in (at present) the following possible ways:
//
// (a) as the identity function x |--> x
//     [constructor f(IDENTITY)]
//
// (b) as the function x |--> x + n, where n is a constant 
//     [constructor f(ADD_CONST, n)]
//
// (c) as the function x |--> v[x], where v is an externally maintained vec<int>
//     [constructor f(APPLY_EXTERNAL_VEC, v)].
//
// To call the function, use f(x).

#ifndef MAP_INT_INT_H
#define MAP_INT_INT_H

#include "CoreTools.h"

const int UNDEFINED = -1;
const int IDENTITY = 0;
const int ADD_CONST = 1;
const int APPLY_EXTERNAL_VEC = 2;

class map_int_int {

     public:

     // ----------------------------------------------------------------------------

     map_int_int( )
     {    fun_type_ = UNDEFINED;    }

     // ----------------------------------------------------------------------------

     void Def( int fun_type )
     {    ForceAssertEq( fun_type, IDENTITY );
          fun_type_ = fun_type;    }

     map_int_int( int fun_type ) { Def(fun_type); }

     // ----------------------------------------------------------------------------

     void Def( int fun_type, int n )
     {    ForceAssertEq( fun_type, ADD_CONST );
          fun_type_ = fun_type;
          to_add_ = n;    }

     map_int_int( int fun_type, int n ) { Def(fun_type, n); }

     // ----------------------------------------------------------------------------

     void Def( int fun_type, const vec<int>& v )
     {    ForceAssertEq( fun_type, APPLY_EXTERNAL_VEC );
          fun_type_ = fun_type;
          vec_to_apply_ = &v;    }

     map_int_int( int fun_type, const vec<int>& v ) { Def(fun_type, v); }

     // ----------------------------------------------------------------------------

     int operator( )( int x ) const
     {    if ( fun_type_ == IDENTITY ) return x;
          else if ( fun_type_ == ADD_CONST ) return x + to_add_;
          else // ( fun_type_ == APPLY_EXTERNAL_VEC )
          {    return (*vec_to_apply_)[x];    }    }

     private:

     int fun_type_;
     int to_add_;
     const vec<int>* vec_to_apply_;

};

class map_longlong_longlong {

     public:

     // ----------------------------------------------------------------------------

     map_longlong_longlong( )
     {    fun_type_ = UNDEFINED;    }

     // ----------------------------------------------------------------------------

     void Def( longlong fun_type )
     {    ForceAssertEq( fun_type, IDENTITY );
          fun_type_ = fun_type;    }

     map_longlong_longlong( longlong fun_type ) { Def(fun_type); }

     // ----------------------------------------------------------------------------

     void Def( longlong fun_type, longlong n )
     {    ForceAssertEq( fun_type, ADD_CONST );
          fun_type_ = fun_type;
          to_add_ = n;    }

     map_longlong_longlong( longlong fun_type, longlong n ) { Def(fun_type, n); }

     // ----------------------------------------------------------------------------

     void Def( longlong fun_type, const vec<longlong>& v )
     {    ForceAssertEq( fun_type, APPLY_EXTERNAL_VEC );
          fun_type_ = fun_type;
          vec_to_apply_ = &v;    }

     map_longlong_longlong( longlong fun_type, const vec<longlong>& v ) { Def(fun_type, v); }

     // ----------------------------------------------------------------------------

     longlong operator( )( longlong x ) const
     {    if ( fun_type_ == IDENTITY ) return x;
          else if ( fun_type_ == ADD_CONST ) return x + to_add_;
          else // ( fun_type_ == APPLY_EXTERNAL_VEC )
          {    return (*vec_to_apply_)[x];    }    }

     private:

     longlong fun_type_;
     longlong to_add_;
     const vec<longlong>* vec_to_apply_;

};

#endif
