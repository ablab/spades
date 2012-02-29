///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef MISC
#define MISC

//  Header: /seq/wga/ArachneCVS/Arachne/Misc.h,v 1.9 2000/07/10 13:21:03 jaffe Exp $             
//  #include "mymins.h"  - moved to math/Functions.h 

// The following Macro will prevent gcc version 2.x from issuing an "unused"
// variable waring when we put "$Id: MainArrays.h,v 1.20 2000/06/26 21:26:20 serafim
// Exp $" into each source file.


#ifndef  __DECCXX_VER
//  This is an ugly way to avoid the gcc warning, but it works.
static int dummyfcn( char *qwer, int what ) { return 1; }
#define USE(var) static int use_##var() { return dummyfcn( var, use_##var() ) ; }
// The below is deprecated and gives a warning in gcc 4.2
// static char* dummyvar = "";
static char dummyvar[] = "";
USE(dummyvar)
#else
#define USE(var)
#endif


enum ChainEnd { END_OF_CHAIN     = 1, 
		BROKEN_CHAIN     = 2, 
		THIS_CONTIG      = 3, 
		DIFFERENT_CONTIG = 4, 
		NOT_SET          = 5} ;
 
enum ContigMark { MARK_CONTIGS = 1, 
		  JUST_LOOKING = 2 } ;

enum direction { LOOK_RIGHT = 1, LOOK_LEFT = 2 } ;
inline direction Opposite( direction dir ) 
{ return (dir==LOOK_RIGHT) ? LOOK_LEFT : LOOK_RIGHT ; }

const int max_read_length = 5000 ;
const int HUGE_NEGATIVE_NUMBER = -1000000 ; 
const int Infinity = 1000000006;

const int masks[32] = { 0,        1,        3,        7,         15,        31,        63,         127,        
			255,      511,      1023,     2047,      4095,      8191,      16383,      32767,
			65535,    131071,   262143,   524287,    1048575,   2097151,   4194303,    8388607,
			16777215, 33554431, 67108863, 134217727, 268435455, 536870911, 1073741823, 2147483647 };


const int NO_NEIGHBOR = - 1; 

#ifndef MYBOOL
#define MYBOOL
typedef unsigned char Bool;

const unsigned char True  = 1;
const unsigned char False = 0;


#endif

char PrintBool( Bool b );
int  BoolToInt( Bool b );

#endif
