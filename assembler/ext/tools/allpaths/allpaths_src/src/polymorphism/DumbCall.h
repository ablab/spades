/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// A dumbcall represents what you get if you lay reads on a reference genome and
// call bases that you see at a particular position - either A, C, G, or T, or
// a deletion (D), or insertion (I).

#ifndef DUMB_CALL_H
#define DUMB_CALL_H

#include "math/Functions.h"
#include "dna/Bases.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"

inline char AsBase( const int x )
{    ForceAssertGe( x, 0 );
     ForceAssertLt( x, 6 );
     if ( x == 0 ) return 'A';
     else if ( x == 1 ) return 'C';
     else if ( x == 2 ) return 'G';
     else if ( x == 3 ) return 'T';
     else if ( x == 4 ) return 'D';
     else return 'I';    }

inline char AsNumber( const char x )
{    
     if      ( x == 'A' ) return 0;
     else if ( x == 'C' ) return 1;
     else if ( x == 'G' ) return 2;
     else if ( x == 'T' ) return 3;
     else if ( x == 'D' ) return 4;
     else if ( x == 'I' ) return 5;
     else FatalErr("Couldn't convert call char to a number.");
     return 0;
}

class dumbcall {

     public:

     int base[6];

     dumbcall( ) { memset(base, 0x00, sizeof(int) * 6); }

     dumbcall( const int A, const int C, const int G, const int T, const int D,
          const int I ) 
     {    base[0] = A;
          base[1] = C;
          base[2] = G;
          base[3] = T;
          base[4] = D;
          base[5] = I;    }

     int A( ) const { return base[0]; }
     int C( ) const { return base[1]; }
     int G( ) const { return base[2]; }
     int T( ) const { return base[3]; }
     int D( ) const { return base[4]; }
     int I( ) const { return base[5]; }

     double BestAltFrac( const char refbase )
     {    if ( CountAll( ) == 0 ) return 0;
          int M = 0;
          for ( int j = 0; j < 6; j++ )
               if ( j != refbase ) M = Max( M, base[j] );
          return double(M)/double( CountAll( ) );    }

     // PolyCI: return True if there is a base other than refbase such that
     // we are confident with probability p that the true fraction of the base
     // is at least minfrac.
         
     Bool PolyCI( const int refbase, const double p, const double minfrac ) const
     {    int total = CountAll( );
          for ( int j = 0; j < 6; j++ )
          {    if ( j == refbase ) continue;
               double z = InverseNormalCDF( 1.0 - (1.0-p)/2.0 );
               double low = base[j] - z * sqrt(base[j]);
               if ( low / double(total) >= minfrac ) return True;    }
          return False;    }
         
     Bool PolyCINoIndels( 
          const int refbase, const double p, const double minfrac ) const
     {    int total = CountAll( );
          for ( int j = 0; j < 4; j++ )
          {    if ( j == refbase ) continue;
               double z = InverseNormalCDF( 1.0 - (1.0-p)/2.0 );
               double low = base[j] - z * sqrt(base[j]);
               if ( low / double(total) >= minfrac ) return True;    }
          return False;    }
         
     Bool Poly( const int refbase, const double frac, const int mincount = 0 ) const
     {    for ( int j = 0; j < 6; j++ )
          {    if ( j == refbase ) continue;
               if ( base[j] < mincount ) continue;
               if ( double(base[j]) > double(base[refbase])*frac ) return True;    }
          return False;    }
         
     Bool PolyNoIndels( const int refbase, const double frac, 
          const int mincount = 0 ) const
     {    for ( int j = 0; j < 4; j++ )
          {    if ( j == refbase ) continue;
               if ( base[j] < mincount ) continue;
               if ( double(base[j]) > double(base[refbase])*frac ) return True;    }
          return False;    }
         
     String PolyNoIndelsCalls( const int refbase, const double frac, 
          const int mincount = 0 ) const
     {    String answer;
          for ( int j = 0; j < 4; j++ )
          {    if ( j == refbase ) continue;
               if ( base[j] < mincount ) continue;
               if ( double(base[j]) > double(base[refbase])*frac ) 
                    answer += as_base(j);    }
          return answer;    }
         
     Bool Poly( const int refbase, const double frac ) const
     {    for ( int j = 0; j < 6; j++ )
          {    if ( j == refbase ) continue;
               if ( double(base[j]) > double(base[refbase])*frac ) return True;    }
          return False;    }

     char PolyBase( const int refbase, const double frac ) const
     {    for ( int j = 0; j < 6; j++ )
          {    if ( j == refbase ) continue;
               if ( double(base[j]) > double(base[refbase])*frac ) return True;    }
          return -1;    }

     int PolyCount( const int refbase, const double frac ) const
     {    for ( int j = 0; j < 6; j++ )
          {    if ( j == refbase ) continue;
               if ( double(base[j]) > double(base[refbase])*frac ) return True;    }
          return -1;    }

     int Count( const int x ) const
     {    return base[x];    }

     int CountAll( ) const
     {    return base[0]+base[1]+base[2]+base[3]+base[4]+base[5];    }

     int CountOther( const int x ) const
     {    int sum = 0;
          for ( int i = 0; i < 6; i++ )
               if ( i != x ) sum += base[i];
          return sum;    }

     String ShowAll( ) const
     {    String s;
          for ( int i = 0; i < 6; i++ )
          {    if ( base[i] > 0 )
                    s += String(AsBase(i)) + "[" + ToString(base[i]) + "]";    }
          return s;    }

     String ShowAll( const char refbase, const int min_to_show = 1 ) const
     {    String s;
          if ( CountAll( ) > 0 )
          {    s += ToString(base[(int)refbase]);
               for ( int i = 0; i < 6; i++ )
               {    if ( base[i] >= min_to_show && i != refbase )
                    {    s += "," + ToString(base[i]) + "["
                              + String(AsBase(i)) + "]";    }    }    }
          return s;    }

     String ShowAll2( char refbase, const int min_to_show = 1 ) const
     {
          refbase = AsNumber(refbase);

          char buf[4096];
          sprintf(buf, "%c[%d] ", AsBase(refbase), base[(int)refbase]);
          for (int i = 0; i < 6; i++)
          {
            if (i == refbase) { continue; }
            if (base[i] < min_to_show) { continue; }
            sprintf(buf, "%s%c[%d]", buf, AsBase(i), base[i]);
          }

          String s = buf;
		  return s;
     }

    String TabularString(char refbase)
    {
        char buf[4096];
        sprintf(buf, "%c %d %d %d %d", 
                        refbase,
                        base[0],
                        base[1],
                        base[2],
                        base[3]);
        String s = buf;
        return s;
    }

     dumbcall& operator+=( const dumbcall& c )
     {    for ( int i = 0; i < 6; i++ )
               this->base[i] += c.base[i];
          return *this;    }
};

TRIVIALLY_SERIALIZABLE(dumbcall);

class dumbercall {

     public:

     int base[4];

     dumbercall( )
     {    for ( int i = 0; i < 6; i++ )
               base[i] = 0;    }

     dumbercall( const int A, const int C, const int G, const int T )
     {    base[0] = A;
          base[1] = C;
          base[2] = G;
          base[3] = T;    }

     int A( ) const { return base[0]; }
     int C( ) const { return base[1]; }
     int G( ) const { return base[2]; }
     int T( ) const { return base[3]; }
         
     Bool Poly( const int refbase, const double frac, const int mincount ) const
     {    for ( int j = 0; j < 4; j++ )
          {    if ( j == refbase ) continue;
               if ( base[j] < mincount ) continue;
               if ( double(base[j]) > double(base[refbase])*frac ) return True;    }
          return False;    }
         
     int Count( const int x ) const
     {    return base[x];    }

     int CountAll( ) const
     {    return base[0]+base[1]+base[2]+base[3];    }

     int CountOther( const int x ) const
     {    int sum = 0;
          for ( int i = 0; i < 4; i++ )
               if ( i != x ) sum += base[i];
          return sum;    }

     String ShowOther( const int x ) const
     {    String s;
          for ( int i = 0; i < 4; i++ )
          {    if ( i != x && base[i] > 0 )
                    s += String(as_base(i)) + "[" + ToString(base[i]) + "]";    }
          return s;    }
               
};

typedef SerfVec<dumbcall> DumbcallVec;
typedef MasterVec<DumbcallVec> VecDumbcallVec;


// Implement A = A + B
void SumDumbCalls(VecDumbcallVec& A, VecDumbcallVec& B);

#endif
