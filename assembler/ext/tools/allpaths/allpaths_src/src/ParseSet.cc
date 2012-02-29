///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "ParseSet.h"
#include "Set.h"
#include "random/Random.h"

#define FLUNK( routine )                                                      \
{    if (ABORT_IF_BAD)                                                        \
          FatalErr( #routine << ": illegal descriptor: " << descrip << "." ); \
     status = 1;                                                              \
     return;    }

void ParseIntSet( String descrip, vec<int>& answer, int& status, 
                  Bool ABORT_IF_BAD, bool sortAnswer, const int start,
                  const int stop ) {
  vec<longlong> lanswer;
  ParseLongLongSet(descrip, lanswer, status, ABORT_IF_BAD, sortAnswer, start, stop);
  for (int i = 0; i < lanswer.isize(); ++i) {
    if ((lanswer[i] <= numeric_limits<int>::max()) && 
	lanswer[i] >= numeric_limits<int>::min())
      answer.push_back(static_cast<int>(lanswer[i]));
    else if (ABORT_IF_BAD)
      FatalErr( "ParseIntSet: descriptor out of range: " << descrip << "." );
  }
  return;
}

void TestStartStop( const vec<longlong>& answer, int& status, Bool ABORT_IF_BAD, 
     const longlong start, const longlong stop )
{    Bool start_stop_provided = ( !( start == 0 && stop == -1 ) );
     if (start_stop_provided)
     {    for ( ulonglong j = 0; j < answer.size( ); j++ )
          {    if ( !( answer[j] >= start && answer[j] < stop ) )
               {    status = 1;
                    if (ABORT_IF_BAD) 
                    {    PRINT3( start, stop, answer[j] );
                         FatalErr( "[start,stop) range violated" );    }
                    return;    }    }    }    }

void ParseLongLongSet( String descrip, vec<longlong>& answer, int& status, 
		       Bool ABORT_IF_BAD, bool sortAnswer, const longlong start,
                       const longlong stop )
{    
     Bool start_stop_provided = ( !( start == 0 && stop == -1 ) );

     if ( descrip.size( ) == 0  || descrip == "{}" )
       return;

     if ( descrip.Contains( "random:", 0 ) && descrip.After( "random:" ).IsInt( ) )
     {    if ( !start_stop_provided )
          {    status = 1;
               if (ABORT_IF_BAD) 
               {    FatalErr( "random: can only be used if start "
                         "and stop are specified." );    }
               return;    }
          longlong r = descrip.After( "random:" ).Int( );
          if ( r > stop - start )
          {    status = 1;
               if (ABORT_IF_BAD) FatalErr( "random:n - n is too large." );
               return;    }
          set<longlong> R;
          for ( longlong i = 0; i < r; i++ )
          {    longlong x = ( big_random( ) % (stop - start) ) + start;
               if ( Member( R, x ) ) i--;
               else R.insert(x);    }
          for ( set<longlong>::iterator i = R.begin( ); i != R.end( ); i++ )
               answer.push_back(*i);
          TestStartStop( answer, status, ABORT_IF_BAD, start, stop );
          status = 0;
          return;    }

     if ( descrip.Contains( "|" ) )
     {    String descrip1 = descrip.Before( "|" ), descrip2 = descrip.After( "|" );
          vec<longlong> answer1, answer2;
          int status1, status2;
          ParseLongLongSet( descrip1, answer1, status1, ABORT_IF_BAD );
          ParseLongLongSet( descrip2, answer2, status2, ABORT_IF_BAD );
          if ( status1 == 1 || status2 == 1 ) status = 1;
          else 
          {    status = 0;
               answer = answer1;
               answer.append(answer2);
               if (sortAnswer) UniqueSort(answer);    }
          TestStartStop( answer, status, ABORT_IF_BAD, start, stop );
          return;    }
     
     if ( descrip.IsInt( ) )
     {    answer.resize(1);
          answer[0] = descrip.Int( );    }

     else if ( descrip.Contains( "{", 0 ) )
     {    if ( !descrip.Contains( "}", -1 ) ) FLUNK(ParseIntSet);
          answer.clear( );
          descrip = descrip.After( "{" );
          while( descrip.Contains( "," ) )
          {    String next = descrip.Before( "," );
               if ( !next.IsInt( ) ) FLUNK(ParseIntSet);
               if ( !next.empty( ) ) answer.push_back( next.Int( ) );
               descrip = descrip.After( "," );    }
          descrip.erase( descrip.size( ) - 1, 1 );
          if ( !descrip.IsInt( ) ) FLUNK(ParseIntSet);
          if ( !descrip.empty( ) ) answer.push_back( descrip.Int( ) );
          if (sortAnswer) UniqueSort(answer);    }

     else if ( descrip.Contains( "[", 0 ) && descrip.Contains( "]", -1 ) )
     {    descrip = descrip.After( "[" );
          if ( !descrip.Contains( "," ) ) FLUNK(ParseIntSet);
          String first = descrip.Before( "," );
          if ( first.empty( ) || !first.IsInt( ) ) FLUNK(ParseIntSet);
          longlong f = first.Int( );
          descrip = descrip.After( "," );
          descrip.erase( descrip.size( ) - 1, 1 );
          if ( descrip.empty( ) || !descrip.IsInt( ) ) FLUNK(ParseIntSet);
          longlong l = descrip.Int( );
          if ( l < f ) FLUNK(ParseIntSet);
          answer.resize( l - f + 1 );
          for ( int i = f; i <= l; i++ )
               answer[ i - f ] = i;    }

     else if ( descrip.Contains( "[", 0 ) && descrip.Contains( ")", -1 ) )
     {    descrip = descrip.After( "[" );
          if ( !descrip.Contains( "," ) ) FLUNK(ParseIntSet);
          String first = descrip.Before( "," );
          if ( first.empty( ) || !first.IsInt( ) ) FLUNK(ParseIntSet);
          longlong f = first.Int( );
          descrip = descrip.After( "," );
          descrip.erase( descrip.size( ) - 1, 1 );
          if ( descrip.empty() || !descrip.IsInt( ) ) FLUNK(ParseIntSet);
          longlong l = descrip.Int( );
          if ( l < f ) FLUNK(ParseIntSet);
          answer.resize( l - f );
          for ( int i = f; i < l; i++ )
               answer[ i - f ] = i;    }

     else if ( descrip.Contains( "@", 0 ) )
     {    descrip = descrip.After( "@" );
          if ( !IsRegularFile(descrip) ) FLUNK(ParseIntSet);
          answer.clear( );
          String entry;
          Ifstream( in, descrip );
          while(1)
          {    in >> entry;
               if ( !in ) break;
               if ( !entry.IsInt( ) ) FLUNK(ParseIntSet);
               answer.push_back( entry.Int( ) );    }
          if (sortAnswer) UniqueSort(answer);    }
               
     else FLUNK(ParseIntSet);

     TestStartStop( answer, status, ABORT_IF_BAD, start, stop );

     status = 0;    }

void ParseIntSet( String descrip, vec<int>& answer, bool sortAnswer,
     const int start, const int stop )
{    int status;
     ParseIntSet( descrip, answer, status, True, sortAnswer, start, stop );    }

void ParseLongLongSet( String descrip, vec<longlong>& answer, bool sortAnswer,
     const longlong start, const longlong stop )
{    int status;
     ParseLongLongSet( descrip, answer, status, True, sortAnswer, start, stop );    }

void ParseStringSet( String descrip, vec<String>& answer, bool recurse )
{    
     answer.clear( );

     // Do case where descrip is empty.

     if ( descrip.size( ) == 0  || descrip == "{}" ) { }

     // Do case where descrip does not contain curly brackets.
       
     else if ( !descrip.Contains( "{" ) ) 
     {
       if ( descrip.Contains( "@", 0 ) )
       {
	 descrip = descrip.After( "@" );
	 if ( !IsRegularFile(descrip) ) 
	 {
	   FatalErr( "ParseStringSet: file not found: " << descrip << "." );
	 } 
	 answer.clear( );
	 String entry;
	 Ifstream( in, descrip );
	 while(1)
          {  
	    in >> entry;
	    if ( !in ) break;
	    answer.push_back( entry );   
	  }
       }
       else
       {
	 answer.push_back(descrip);
       }
     }

     // Do case where descrip does not contain asterisk.

     else if ( !descrip.Contains( "*" ) )
     {    String x = descrip.Before( "{" ), rest = descrip.After( "{" );
          if ( !rest.Contains( "}" ) ) PRINT2( descrip, rest );
          ForceAssert( rest.Contains( "}" ) );
          String middle = rest.RevBefore( "}" ), y = rest.RevAfter( "}" );
          while( middle.nonempty( ) ) 
          {    int comma, brackcount = 0;
               for ( comma = 0; comma < middle.isize( ); comma++ )
               {    if ( middle[comma] == ',' && brackcount == 0 ) break;
                    if ( middle[comma] == '{' ) brackcount++;
                    if ( middle[comma] == '}' ) brackcount--;    }
               if ( comma < middle.isize( ) )
               {    answer.push_back( x + middle.substr( 0, comma ) + y );
                    middle = middle.substr( comma+1, -1 );    }
               else
               {    answer.push_back( x + middle + y );
                    middle = "";    }    }    }

     // Do the other cases.

     else 
     {    String descripx = descrip;
          String descripy;
          int mult1 = 1;
          if ( descrip.Contains( "}*" ) ) mult1 = descrip.After( "}*" ).Int( );
          else ForceAssert( descrip[ descrip.size( ) - 1 ] == '}' );
          descripx = descripx.After( "{" ).RevBefore( "}" );
          while( descripx.nonempty( ) ) 
          {    int comma, brackcount = 0;
               for ( comma = 0; comma < descripx.isize( ); comma++ )
               {    if ( descripx[comma] == ',' && brackcount == 0 ) break;
                    if ( descripx[comma] == '{' ) brackcount++;
                    if ( descripx[comma] == '}' ) brackcount--;    }
               if ( comma < descripx.isize( ) )
               {    descripy = descripx.substr( 0, comma );
                    descripx = descripx.substr( comma + 1, -1 );    }
               else 
               {    descripy = descripx;
	            descripx = "";    }
               int mult2 = 1;
               if ( descripy.Contains("*") && descripy.After("*").IsInt() ) 
               {    mult2 = descripy.After("*").Int();
	            descripy = descripy.Before("*");    }
               for ( int m = 0; m < mult1 * mult2; m++ )
	            answer.push_back( descripy );    }    }    

     if ( recurse ) {
       // if we did anything at this level, recurse again
       if ( answer.size() > 1 || answer.front() != descrip ) {
         vec<String> recursive_answer;
         for ( unsigned int i = 0; i < answer.size(); ++i ) {
           vec<String> one_answer;
           ParseStringSet( answer[i], one_answer, recurse );
           recursive_answer.append( one_answer );
         }
         answer.swap( recursive_answer );
       }
     }
}

void ParseDoubleSet( const String & descrip, vec<double> & answer, bool sortAnswer ) {
  Bool ABORT_IF_BAD = True;
  int status = -1;

  answer.clear( );

  // Do case where descrip is empty.

  if ( descrip.size( ) == 0  || descrip == "{}" )
    return;

  if ( descrip[0] == '@' ) {
    String filename = descrip.After( "@" );
    if ( !IsRegularFile(filename) ) FLUNK(ParseDoubleSet);
    answer.clear( );
    String entry;
    Ifstream( in, filename );
    while(1) {
      in >> entry;
      if ( !in ) break;
      if ( !entry.IsDouble( ) ) FLUNK(ParseDoubleSet);
      answer.push_back( entry.Double( ) ); 
    }
    if (sortAnswer) UniqueSort(answer);    
  }
  else if ( descrip.IsDouble() ) {
    answer.push_back( descrip.Double( ) );
  }
  else {
    String descripx = descrip;
    descripx = descripx.SafeAfter( "{" ).SafeBefore( "}" );
    while(1) {    
      if ( descripx.Contains( "," ) ) {    
        answer.push_back( descripx.Before( "," ).Double() );
        descripx = descripx.After( "," );    
      }
      else {    
        answer.push_back(descripx.Double());
        break;    
      }    
    }

    if (sortAnswer) {
      sort(answer.begin(), answer.end());
      answer.erase(unique(answer.begin(), answer.end()), answer.end());
    }
  }
}   
  
