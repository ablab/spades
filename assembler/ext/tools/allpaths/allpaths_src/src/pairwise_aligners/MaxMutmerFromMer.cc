// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#include "Basevector.h"
#include "pairwise_aligners/MaxMutmerFromMer.h"

// ================================================================================
//
// MaxMutmerFromMer: Given two reads (rd1, rd2), and positions on both of them
// (pos1, pos2), and a match length (len), corresponding to a "len-mer", extend
// on both ends to obtain a max-mutmer.
//
// This is a time-critical piece of code.
//
// MaxMutmerFromMerRev: same, but use rc(rd2) instead of rd2.
//
// ================================================================================

void MaxMutmerFromMer( int& pos1, int& pos2, int& len, int& errors,
     const basevector& rd1, const basevector& rd2, Bool strict )
{    register unsigned int xpos1 = pos1, xpos2 = pos2, xlen = len, xerrors = 0;

     // Extend on left side.

     int goodlength = xlen; // Any big number.
     while(1)
     {    if ( xpos1 == 0 || xpos2 == 0 ) break;
          if( rd1[xpos1-1] == rd2[xpos2-1] )
          {    --xpos1; --xpos2;
               ++xlen;    
               ++goodlength;    }
          else if (strict) break;
          else if ( goodlength <= 1 ) break;
          else
          {    --xpos1; --xpos2;
               ++xlen;
               ++xerrors;
               goodlength = 0;    }    }
      if ( goodlength <= 1 ) // Back up.
      {    xpos1 += goodlength + 1; xpos2 += goodlength + 1;
           xlen -= goodlength + 1;    
           --xerrors;    }

     // Extend on right side.

     goodlength = xlen;
     while(1)
     {    if ( xpos1+xlen >= rd1.size( ) || xpos2+xlen >= rd2.size( ) ) break;
          if ( rd1[xpos1+xlen] == rd2[xpos2+xlen] )
          {    ++xlen;
               ++goodlength;    }
          else if (strict) break;
          else if ( goodlength <= 1 ) break;
          else
          {    ++xlen;
               ++xerrors;
               goodlength = 0;    }    }
     if ( goodlength <= 1 )
     {    xlen -= goodlength + 1;
          --xerrors;    }

     if ( xerrors > 31 ) xerrors = 31;

     pos1 = xpos1; pos2 = xpos2;
     len = xlen;
     errors = xerrors;    }

void MaxMutmerFromMerRev( int& pos1, int& pos2, int& len, int& errors,
     const basevector& rd1, const basevector& rd2, Bool strict )
{    register unsigned int xpos1 = pos1, xpos2 = pos2, xlen = len, xerrors = 0;

     // Extend on left side.

     int goodlength = xlen; // Any big number.
     while(1)
     {    if ( xpos1 == 0 || xpos2 == 0 ) break;
          if( rd1[xpos1-1] == 3 - rd2[rd2.size( ) - xpos2] )
          {    --xpos1; --xpos2;
               ++xlen;    
               ++goodlength;    }
          else if (strict) break;
          else if ( goodlength <= 1 ) break;
          else
          {    --xpos1; --xpos2;
               ++xlen;
               ++xerrors;
               goodlength = 0;    }    }
      if ( goodlength <= 1 ) // Back up.
      {    xpos1 += goodlength + 1; xpos2 += goodlength + 1;
           xlen -= goodlength + 1;    
           --xerrors;    }

     // Extend on right side.

     goodlength = xlen;
     while(1)
     {    if ( xpos1+xlen >= rd1.size( ) || xpos2+xlen >= rd2.size( ) ) break;
          if ( rd1[xpos1+xlen] == 3 - rd2[rd2.size( ) - 1 - (xpos2+xlen) ] )
          {    ++xlen;
               ++goodlength;    }
          else if (strict) break;
          else if ( goodlength <= 1 ) break;
          else
          {    ++xlen;
               ++xerrors;
               goodlength = 0;    }    }
     if ( goodlength <= 1 )
     {    xlen -= goodlength + 1;
          --xerrors;    }

     if ( xerrors > 31 ) xerrors = 31;

     pos1 = xpos1; pos2 = xpos2;
     len = xlen;
     errors = xerrors;    }
