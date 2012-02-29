///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "CoreTools.h"
#include "Fastavector.h"
#include "PackAlign.h"
#include "ShortVector.h"
#include "VecTemplate.h"
#include "feudal/OuterVecDefs.h"
#include "feudal/SmallVecDefs.h"
#include "math/Functions.h"

int PosDel( int x )
{
  return (x < 0) ? x : (x - 1);    
}

int PosDel2( int x )
{
  return PosDel(x) & Bits2;    
}

int PosDel4( int x )
{
  return PosDel(x) & Bits4;    
}

int PosAdd( int x )
{
  return (x < 0) ? x : (x + 1);    
}

int PosAdd2( int x )
{
  if ( x <= 1 ) return x + 1;
  return x - 4;    
}

void packalign::ConstructorCore( int pos1, int pos2, 
                                 const avector<int>& gaps, const avector<int>& lengths, 
                                 int nblocks )
{
    
  int n = gaps.length;
  if ( nblocks >= 0 ) n = nblocks;

#ifndef NDEBUG
  {
    for ( int i = 0; i < n; i++ )
      Assert( lengths(i) >= 0 );    
  }
#endif

  // gaps (after the 0th) shouldn't be zero, but just in case, we don't
  // want to create garbage, so we force type 3

  Bool zero_gap = False;
  for ( int i = 1; i < n; i++ )
    if ( gaps(i) == 0 )
    {
      zero_gap = True;
      break;    
    }

  if ( pos1 <= 1023 && pos2 <= 1023 && n <= 6 && n >= 1 
       && gaps(0) == 0 && pos1 >= 0 && pos2 >= 0 && !zero_gap )
  {
    int i;
    for ( i = 0; i < n; i++ )
      if ( lengths(i) > 1023 || Abs(gaps(i)) > 2 ) break;
    if ( i == n )
    {
      ConvertToType0( pos1, pos2, gaps, lengths, nblocks );
      Assert( Control( ) == 0 ); // XXX
      return;    
    }   
  }

  if ( pos1 <= 4095 && pos2 <= 4095 && n >= 1 && gaps(0) == 0 
       && pos1 >= 0 && pos2 >= 0 && !zero_gap )
  {
    int i;
    for ( i = 0; i < n; i++ )
      if ( lengths(i) > 4095 || Abs(gaps(i)) > 8 ) break;
    if ( i == n )
    {
      ConvertToType1( pos1, pos2, gaps, lengths, nblocks );
      Assert( Control( ) == 1 ); // XXX
      return;    
    }   
  }

  ConvertToType2( pos1, pos2, gaps, lengths, nblocks );    
  Assert( Control( ) == 2 ); // XXX

}

packalign::packalign( int pos1, int pos2, 
                      const avector<int>& gaps, const avector<int>& lengths )
{
  ForceAssert( gaps.x != 0 );
  ForceAssert( lengths.x != 0 );
  ForceAssert( gaps.length <= 65535 );
  ForceAssert( gaps.length == lengths.length );
  ConstructorCore( pos1, pos2, gaps, lengths );    
  /* // XXX
     int actual_pos1, actual_pos2; // XXX
     avector<int> actual_gaps, actual_lengths; // XXX
     Unpack( actual_pos1, actual_pos2, actual_gaps, actual_lengths ); // XXX
     for ( unsigned int i = 0; i < gaps.length; i++ ) // XXX
     if ( gaps(i) != actual_gaps(i) ) // XXX
     {
     cout << "packalign constructor failure\n"; // XXX
     PRINT3( i, gaps(i), actual_gaps(i) ); // XXX
     exit(1);    
     } // XXX
     for ( unsigned int i = 0; i < lengths.length; i++ ) // XXX
     if ( lengths(i) != actual_lengths(i) ) // XXX
     {
     cout << "packalign constructor failure\n"; // XXX
     PRINT3( i, lengths(i), actual_lengths(i) ); // XXX
     exit(1);    
     } // XXX
     if ( pos1 != actual_pos1 || pos2 != actual_pos2 ) // XXX
     {
     cout << "packalign constructor failure\n"; // XXX
     PRINT2( pos1, actual_pos1 ); // XXX
     PRINT2( pos2, actual_pos2 ); // XXX
     exit(1);    
     } // XXX
  */ // XXX

}

packalign::packalign( const align& a )
{
  ConstructorCore( a.pos1( ), a.pos2( ), a.Gaps( ), a.Lengths( ), 
                   a.Nblocks( ) );    
}

packalign::packalign( const packalign& p )
{
    
  word_[0] = p.word_[0];
  int control = p.Control( );
  if ( control == 0 )
  {
    word_[1] = p.word_[1];
    word_[2] = p.word_[2];    
  }
  else if ( control == 1 )
  {
    short_pointer_or_words pw, pw_new;
    pw.x[0] = p.word_[1];
    pw.x[1] = p.word_[2];
    unsigned short* ptr = pw.p;
    int nwords = 1 + ptr[0];
    unsigned short* ptr_new = new unsigned short[nwords];
    for ( int i = 0; i < nwords; i++ )
      ptr_new[i] = ptr[i];
    pw_new.p = ptr_new;
#ifdef IS_32_BIT
    word_[1] = (unsigned int) ptr_new;
#else 
    word_[1] = pw_new.x[0];
    word_[2] = pw_new.x[1];
#endif
  }
  else if ( control == 2 )
  {
    int_pointer_or_words pw, pw_new;
    pw.x[0] = p.word_[1];
    pw.x[1] = p.word_[2];
    unsigned int* ptr = pw.p;
    int nwords = 2 + 2 * (p.word_[0] & Bits16);
    unsigned int* ptr_new = new unsigned int[nwords];
    for ( int i = 0; i < nwords; i++ )
      ptr_new[i] = ptr[i];
    pw_new.p = ptr_new;
#ifdef IS_32_BIT
    word_[1] = (unsigned int) ptr_new;
#else 
    word_[1] = pw_new.x[0];
    word_[2] = pw_new.x[1];
#endif
  }
  else if ( control == 7 ) word_[0] = 7u << 29;    
}

// int pos1, pos2;
// avector<int> gaps, lengths;
// p.Unpack( pos1, pos2, gaps, lengths );
// ConstructorCore( pos1, pos2, gaps, lengths );    }

void packalign::Set( int pos1, int pos2, 
                     const avector<int>& gaps, const avector<int>& lengths, int nblocks )
{
  int n = gaps.length;
  if ( nblocks >= 0 ) n = nblocks;
  else ForceAssert( gaps.length == lengths.length );

  // This actually asserts!!!
  // if ( n > 1000 ) // XXX
  // {    PRINT2(pos1, pos2); // XXX
  //      PRINT(n); // XXX
  //      Assert( 0 == 1 );    }
  // XXX

  ForceAssert( gaps.x != 0 );
  ForceAssert( lengths.x != 0 );
  ForceAssert( n <= 65535 );

  int control = Control( );
  if ( control == 0 )      DeleteType0( );
  else if ( control == 1 ) DeleteType1( );
  else if ( control == 2 ) DeleteType2( );

  ConstructorCore( pos1, pos2, gaps, lengths, n );    
}

packalign::~packalign( )
{
  int control = Control( );
  if ( control == 0 )      DeleteType0( );
  else if ( control == 1 ) DeleteType1( );
  else if ( control == 2 ) DeleteType2( );    
}

packalign& packalign::operator=( const packalign& p )
{
  // Check for the dangerous case where *this == p, when control = 1 or 2. 
  // Also deal with the harmless case of equality in case control = 0.

  if ( word_[0] == p.word_[0] && word_[1] == p.word_[1] 
       && word_[2] == p.word_[2] ) return *this;

  // Here's the simple version.  The fancy version is only a little bit faster.
  //
  //      int pos1, pos2;
  //      avector<int> gaps, lengths;
  //      p.Unpack( pos1, pos2, gaps, lengths );
  //      Set( pos1, pos2, gaps, lengths );
  //      return *this;

  int control = Control( );
  if ( control == 0 )      DeleteType0( );
  else if ( control == 1 ) DeleteType1( );
  else if ( control == 2 ) DeleteType2( );

  word_[0] = p.word_[0];
  control = p.Control( );
  if ( control == 0 )
  {
    word_[1] = p.word_[1];
    word_[2] = p.word_[2];    
  }
  else if ( control == 1 )
  {
    short_pointer_or_words pw, pw_new;
    pw.x[0] = p.word_[1];
    pw.x[1] = p.word_[2];
    unsigned short* ptr = pw.p;
    int nwords = 1 + ptr[0];
    unsigned short* ptr_new = new unsigned short[nwords];
    for ( int i = 0; i < nwords; i++ )
      ptr_new[i] = ptr[i];
    pw_new.p = ptr_new;
#ifdef IS_32_BIT
    word_[1] = (unsigned int) ptr_new;
#else 
    word_[1] = pw_new.x[0];
    word_[2] = pw_new.x[1];
#endif
  }
  else if ( control == 2 )
  {
    int_pointer_or_words pw, pw_new;
    pw.x[0] = p.word_[1];
    pw.x[1] = p.word_[2];
    unsigned int* ptr = pw.p;
    int nwords = 2 + 2 * (p.word_[0] & Bits16);
    unsigned int* ptr_new = new unsigned int[nwords];
    for ( int i = 0; i < nwords; i++ )
      ptr_new[i] = ptr[i];
    pw_new.p = ptr_new;
#ifdef IS_32_BIT
    word_[1] = (unsigned int) ptr_new;
#else
    word_[1] = pw_new.x[0];
    word_[2] = pw_new.x[1];
#endif
  }

  return *this;    
}

void packalign::ConvertToType0( int pos1, int pos2, 
                                const avector<int>& gaps, const avector<int>& lengths, 
                                int nblocks )
{
  int n = gaps.length;

  if ( nblocks >= 0 ) n = nblocks;

  word_[0] = (pos1 << (10 + 3 + 2 + 2 + 2)) ^ (pos2 << (3 + 2 + 2 + 2)) ^
    (n << (2 + 2 + 2));
  if ( n > 1 ) word_[0] ^= ( PosDel2(gaps(1)) << (2 + 2) );
  if ( n > 2 ) word_[0] ^= ( PosDel2(gaps(2)) << 2 );
  if ( n > 3 ) word_[0] ^= PosDel2(gaps(3));
  word_[1] = ( lengths(0) << (10 + 10) );
  if ( n > 4 ) word_[1] ^= ( PosDel2(gaps(4)) << (10 + 10 + 10) );
  if ( n > 1 ) word_[1] ^= ( lengths(1) << 10 );
  if ( n > 2 ) word_[1] ^= lengths(2);
  word_[2] = 0;
  if ( n > 5 ) word_[2] ^= ( PosDel2(gaps(5)) << (10 + 10 + 10) );
  if ( n > 3 ) word_[2] ^= ( lengths(3) << (10 + 10) );
  if ( n > 4 ) word_[2] ^= ( lengths(4) << 10 );
  if ( n > 5 ) word_[2] ^= lengths(5);    
}

void packalign::Unpack( int& pos1, int& pos2, avector<int>& gaps, 
                        avector<int>& lengths ) const
{
  int control = Control( );
  int n = -1;
  Assert( control == 0 || control == 1 || control == 2 || control == 7 );
  if ( control == 7 ) 
  {
    pos1 = 0;
    pos2 = 0;
    gaps.Setsize(0);
    lengths.Setsize(0);    
  }
  else if ( control == 0 ) Unpack0( pos1, pos2, gaps, lengths, n );
  else if ( control == 1 ) Unpack1( pos1, pos2, gaps, lengths, n );
  else if ( control == 2 ) Unpack2( pos1, pos2, gaps, lengths, n );    
}

void packalign::Unpack( int& pos1, int& pos2, avector<int>& gaps, 
                        avector<int>& lengths, int& nblocks ) const
{
  int control = Control( );
  nblocks = 0;
  Assert( control == 0 || control == 1 || control == 2 || control == 7 );
  if ( control == 7 ) 
  {
    pos1 = 0;
    pos2 = 0;
    nblocks = 0;    
  }
  else if ( control == 0 ) Unpack0( pos1, pos2, gaps, lengths, nblocks );
  else if ( control == 1 ) Unpack1( pos1, pos2, gaps, lengths, nblocks );
  else if ( control == 2 ) Unpack2( pos1, pos2, gaps, lengths, nblocks );    
}

void packalign::Unpack0( int& pos1, int& pos2, avector<int>& gaps, 
                         avector<int>& lengths, int& n ) const
{
    
  pos1 = word_[0] >> (10 + 3 + 2 + 2 + 2);
  pos2 = (word_[0] >> (3 + 2 + 2 + 2)) & Bits10;
  int nblocks = (word_[0] >> (2 + 2 + 2)) & Bits3;
  if ( n < 0 )
  {
    gaps.Setsize(nblocks);
    lengths.Setsize(nblocks);    
  }
  else
  {
    if ( gaps.x == 0 || nblocks > (int) gaps.length )
    {
      gaps.Setsize(nblocks);
      lengths.Setsize(nblocks);    
    }
    n = nblocks;    
  }
  gaps(0) = 0;
  lengths(0) = (word_[1] >> (10 + 10)) & Bits10;
  if ( nblocks > 1 )
  {
    gaps(1) = PosAdd2((word_[0] >> (2 + 2)) & Bits2);
    lengths(1) = (word_[1] >> 10) & Bits10;    
  }
  if ( nblocks > 2 )
  {
    gaps(2) = PosAdd2((word_[0] >> 2) & Bits2);
    lengths(2) = word_[1] & Bits10;    
  }
  if ( nblocks > 3 )
  {
    gaps(3) = PosAdd2(word_[0] & Bits2);
    lengths(3) = (word_[2] >> (10 + 10)) & Bits10;    
  }
  if ( nblocks > 4 )
  {
    gaps(4) = PosAdd2((word_[1] >> (10 + 10 + 10)) & Bits2);
    lengths(4) = (word_[2] >> 10) & Bits10;    
  }
  if ( nblocks > 5 )
  {
    gaps(5) = PosAdd2((word_[2] >> (10 + 10 + 10)) & Bits2);
    lengths(5) = word_[2] & Bits10;    
  }    
}

void packalign::ConvertToType1( int pos1, int pos2, 
                                const avector<int>& gaps, const avector<int>& lengths, 
                                int nblocks )
{
  int n = gaps.length;
  if ( nblocks >= 0 ) n = nblocks;
  word_[0] = (1 << 29) ^ (pos1 << 12) ^ pos2;
  unsigned short* ptr = new unsigned short[1 + 2 * n];
  short_pointer_or_words pw;
  pw.p = ptr;

  // If ptr is 4 bytes, put it in word_[1].  Otherwise, put it in word_[1]
  // and word_[2].

#ifdef IS_32_BIT
  word_[1] = (unsigned int) ptr;
#else 
  word_[1] = pw.x[0];
  word_[2] = pw.x[1];    
#endif

  ptr[0] = n;
  for ( int i = 0; i < n; i++ )
    ptr[ 1 + i ] = (PosDel4(gaps(i)) << 12) + lengths(i);    
}

void packalign::Unpack1( int& pos1, int& pos2, avector<int>& gaps, 
                         avector<int>& lengths, int& n ) const
{
  pos1 = (word_[0] >> 12) & Bits12;
  pos2 = word_[0] & Bits12;
  short_pointer_or_words pw;
  pw.x[0] = word_[1];
  pw.x[1] = word_[2];
  unsigned short* ptr = pw.p;

  int nblocks = ptr[0];

  if ( n < 0 )
  {
    gaps.Setsize(nblocks);
    lengths.Setsize(nblocks);    
  }
  else
  {
    if ( gaps.x == 0 || nblocks > (int) gaps.length )
    {
      gaps.Setsize(nblocks);
      lengths.Setsize(nblocks);    
    }
    n = nblocks;    
  }

  for ( int i = 0; i < nblocks; i++ )
  {
    if ( i == 0 ) gaps(i) = 0;
    else gaps(i) = PosAdd( ((short) ptr[1 + i]) >> 12);
    lengths(i) = ptr[1 + i] & Bits12;    
  }    
}

void packalign::ConvertToType2( int pos1, int pos2, 
                                const avector<int>& gaps, const avector<int>& lengths, 
                                int nblocks )
{
  int n = gaps.length;
  if ( nblocks >= 0 ) n = nblocks;
  word_[0] = (2 << 29) ^ n;
  unsigned int* ptr = new unsigned int[2 + 2 * n];
  int_pointer_or_words pw;
  pw.p = ptr;
#ifdef IS_32_BIT
  word_[1] = (unsigned int) ptr;
#else 
  word_[1] = pw.x[0];
  word_[2] = pw.x[1];    
#endif
  ptr[0] = pos1;
  ptr[1] = pos2;
  for ( int i = 0; i < n; i++ )
  {
    ptr[ 2 + 2*i ] = lengths(i);
    ptr[ 3 + 2*i ] = gaps(i);    
  }    
}

void packalign::Unpack2( int& pos1, int& pos2, avector<int>& gaps, 
                         avector<int>& lengths, int& n ) const
{
  int_pointer_or_words pw;
  pw.x[0] = word_[1];
  pw.x[1] = word_[2];
  unsigned int* ptr = pw.p;
  pos1 = ptr[0];
  pos2 = ptr[1];

  int nblocks = word_[0] & Bits16;
  if ( n < 0 )
  {
    gaps.Setsize(nblocks);
    lengths.Setsize(nblocks);    
  }
  else
  {
    if ( gaps.x == 0 || nblocks > (int) gaps.length )
    {
      gaps.Setsize(nblocks);
      lengths.Setsize(nblocks);    
    }
    n = nblocks;    
  }

  for ( int i = 0; i < nblocks; i++ )
  {
    lengths(i) = ptr[ 2 + 2*i ];
    gaps(i) = ptr[ 3 + 2*i ];    
  }    
}

Bool packalign::Dead( ) const
{
  int control = Control( );
  if ( control == 0 ) return ( (word_[0] >> (2 + 2 + 2)) & Bits3 ) == 0;
  if ( control == 1 )
  {
    short_pointer_or_words pw;
    pw.x[0] = word_[1];
    pw.x[1] = word_[2];
    return *pw.p == 0;    
  }
  if ( control == 2 ) return (word_[0] & Bits16) == 0;
  return True;    
}

int packalign::pos1( ) const
{
  align a;
  a.UnpackFrom(*this);
  return a.pos1( );    
}

int packalign::pos2( ) const
{
  align a;
  a.UnpackFrom(*this);
  return a.pos2( );    
}

int packalign::Pos1( ) const
{
  align a;
  a.UnpackFrom(*this);
  return a.Pos1( );    
}

int packalign::Pos2( ) const
{
  align a;
  a.UnpackFrom(*this);
  return a.Pos2( );    
}

void packalign::Kill( )
{
  avector<int> gaps(0), lengths(0);
  Set( 0, 0, gaps, lengths );    
}

void packalign::Flip( )
{
  int pos1, pos2;
  avector<int> gaps, lengths;
  Unpack( pos1, pos2, gaps, lengths );
  for ( unsigned int i = 0; i < gaps.length; i++ )
    gaps(i) = -gaps(i);
  Set( pos2, pos1, gaps, lengths );    
}

packalign packalign::Reverse( int b1_len, int b2_len ) const
{
  int pos1, pos2;
  avector<int> gaps, gaps_new;
  avector<int> lengths, lengths_new;
  Unpack( pos1, pos2, gaps, lengths );
  int n = lengths.length;
  for ( int i = n-1; i >= 0; i-- )
  {
    if ( i == n-1 ) 
    {
      gaps_new.Append(0);
      lengths_new.Append( lengths(n-1) );    
    }
    else 
    {
      gaps_new.Append( gaps(i+1) );
      lengths_new.Append( lengths(i) );    
    }    
  }
  if ( n > 0 && gaps(0) != 0 ) 
  {
    gaps_new.Append( gaps(0) );
    lengths_new.Append(0);    
  }
  int Pos1 = pos1, Pos2 = pos2;
  for ( int j = 0; j < n; j++ )
  {
    if ( gaps(j) > 0 ) Pos2 += gaps(j);
    else if ( gaps(j) < 0 ) Pos1 -= gaps(j);
    Pos1 += lengths(j); 
    Pos2 += lengths(j);    
  }
  return packalign( b1_len - Pos1, b2_len - Pos2, gaps_new, lengths_new );    
}

void packalign::ReverseThis( int b1_len, int b2_len )
{
  *this = Reverse( b1_len, b2_len );    
}

void packalign::SetToFlipOf( const packalign& p )
{
  int pos1, pos2;
  avector<int> gaps, lengths;
  p.Unpack( pos1, pos2, gaps, lengths );
  for ( unsigned int i = 0; i < gaps.length; i++ )
    gaps(i) = -gaps(i);
  Set( pos2, pos1, gaps, lengths );    
}

void packalign::SetToFlipOf( align a )
{
  for ( int i = 0; i < a.Nblocks( ); i++ )
    a.SetGap( i, -a.Gaps(i) );
  Set( a.pos2( ), a.pos1( ), a.Gaps( ), a.Lengths( ), a.Nblocks( ) );    
}

void packalign::SetToReverseFlipOf( const packalign& p, int b1_len, int b2_len )
{
  int pos1, pos2;
  avector<int> gaps, lengths;
  p.Unpack( pos2, pos1, gaps, lengths );
  for ( unsigned int i = 0; i < gaps.length; i++ )
    gaps(i) = -gaps(i);
  int n = lengths.length;
  int new_length = n, ptr = 0;
  if ( n > 0 && gaps(0) != 0 ) ++new_length;
  avector<int> gaps_new(new_length), lengths_new(new_length);
  for ( int i = n-1; i >= 0; i-- )
  {
    if ( i == n-1 ) 
    {
      gaps_new(ptr) = 0;
      lengths_new(ptr) = lengths(n-1);    
    }
    else 
    {
      gaps_new(ptr) = gaps(i+1);
      lengths_new(ptr) = lengths(i);    
    }
    ++ptr;    
  }
  if ( n > 0 && gaps(0) != 0 ) 
  {
    gaps_new(ptr) = gaps(0);
    lengths_new(ptr) = 0;    
  }
  int Pos1 = pos1, Pos2 = pos2;
  for ( int j = 0; j < n; j++ )
  {
    if ( gaps(j) > 0 ) Pos2 += gaps(j);
    else if ( gaps(j) < 0 ) Pos1 -= gaps(j);
    Pos1 += lengths(j); 
    Pos2 += lengths(j);    
  }
  Set( b1_len - Pos1, b2_len - Pos2, gaps_new, lengths_new );    
}

void packalign::SetToReverseFlipOf( align a, int b1_len, int b2_len )
{
  int pos1 = a.pos2( ), pos2 = a.pos1( );
  for ( int i = 0; i < a.Nblocks( ); i++ )
    a.SetGap( i, -a.Gaps(i) );
  int n = a.Nblocks( );
  int new_length = n, ptr = 0;
  if ( n > 0 && a.Gaps(0) != 0 ) ++new_length;
  avector<int> gaps_new(new_length), lengths_new(new_length);
  for ( int i = n-1; i >= 0; i-- )
  {
    if ( i == n-1 ) 
    {
      gaps_new(ptr) = 0;
      lengths_new(ptr) = a.Lengths(n-1);    
    }
    else 
    {
      gaps_new(ptr) = a.Gaps(i+1);
      lengths_new(ptr) = a.Lengths(i);    
    }
    ++ptr;    
  }
  if ( n > 0 && a.Gaps(0) != 0 ) 
  {
    gaps_new(ptr) = a.Gaps(0);
    lengths_new(ptr) = 0;    
  }
  int Pos1 = pos1, Pos2 = pos2;
  for ( int j = 0; j < n; j++ )
  {
    if ( a.Gaps(j) > 0 ) Pos2 += a.Gaps(j);
    else if ( a.Gaps(j) < 0 ) Pos1 -= a.Gaps(j);
    Pos1 += a.Lengths(j); 
    Pos2 += a.Lengths(j);    
  }
  Set( b1_len - Pos1, b2_len - Pos2, gaps_new, lengths_new );    
}

Bool operator==( const packalign& x, const packalign& y )
{
  int xpos1, xpos2, ypos1, ypos2;
  avector<int> xgaps, ygaps; 
  avector<int> xlengths, ylengths;
  x.Unpack( xpos1, xpos2, xgaps, xlengths );
  y.Unpack( ypos1, ypos2, ygaps, ylengths );
  if ( xpos1 != ypos1 || xpos2 != ypos2 ) return False;
  if ( xlengths.length != ylengths.length ) return False;
  for ( int i = 0; i < (int) xlengths.length; i++ )
  {
    if ( xlengths(i) != ylengths(i) ) return False;
    if ( i >= 1 && xgaps(i) != ygaps(i) ) return False;    
  }
  return True;    
}

void align::Compactify( int len1, int len2 )
{
    
  int ni;

  // Clean up alignments that start with gaps.  This should be very rare.
  // The goal here is just to get a valid alignment, not necessarily a good one.

  int first_gap;

  int add_gap_plus = 0;
  int add_gap_minus = 0;

check_first_gap:
  ni = 0;
  if ( nblocks_ == 0 ) goto finish;
  first_gap = gaps_(0);
  if ( first_gap < 0 )
  {
    first_gap = -first_gap; // but on second read
    if ( first_gap <= (int) pos2_ )
    {
      pos2_ -= first_gap;
      lengths_(0) += first_gap;
      gaps_(0) = 0;    
    }
    else
    {
      first_gap -= pos2_;
      gaps_(0) = -first_gap;
      if ( pos2_ > 0 )
      {
        lengths_.Prepend(pos2_);
	gaps_.Prepend(0);
	++nblocks_;    
      }
      pos2_ = 0;    
    }    
  }
  else if ( first_gap > 0 )
  {
    if ( first_gap <= (int) pos1_ )
    {
      pos1_ -= first_gap;
      lengths_(0) += first_gap;
      gaps_(0) = 0;    
    }
    else
    {
      first_gap -= pos1_;
      gaps_(0) = first_gap;
      if ( pos1_ > 0 )
      {
        lengths_.Prepend(pos1_);
	gaps_.Prepend(0);
	++nblocks_;    
      }
      pos1_ = 0;    
    }    
  }
  
  // Remove an initial length of zero, if it exists.
  
  if ( lengths_(0) == 0 ) {    
    if (gaps_(0) > 0)
      add_gap_plus += gaps_(0);
    if (gaps_(0) < 0)
      add_gap_minus += -gaps_(0);

    for ( int i = 1; i < nblocks_; i++ ) { 
      gaps_(i-1) = gaps_(i);
      lengths_(i-1) = lengths_(i);    
    }    
    --nblocks_;
    goto check_first_gap;    
  }
  
  for ( int i = 1; i < nblocks_; i++ )
  {
    if ( lengths_(ni) != 0 && gaps_(i) != 0 )
    {
      ++ni;
      gaps_(ni) = gaps_(i);
      lengths_(ni) = lengths_(i);    
    }
    else
    {
      lengths_(ni) += lengths_(i);
      if ( gaps_(ni) * gaps_(i) < 0 )
	lengths_(ni) += Min( Abs(gaps_(ni)), Abs(gaps_(i)) );
      gaps_(ni) += gaps_(i);    
    }    
  }
  
  if ( lengths_(ni) == 0 ) 
  {
    if ( ni == 0 )
    {
      nblocks_ = 0;
      return;    
    }
    if ( gaps_(ni) < 0 )
    {
      int Pos2 = pos2_;
      for ( int j = 0; j < ni; j++ )
      {
        if ( gaps_(j) > 0 ) Pos2 += gaps_(j);
	Pos2 += lengths_(j);    
      }
      int npos2 = len2 - Pos2;
      Assert( npos2 >= 0 ); // XXX
      if ( npos2 == 0 ) --ni;
      else if ( npos2 >= -gaps_(ni) )
      {
        lengths_(ni-1) -= gaps_(ni);
	--ni;    
      }
      else
      {
        gaps_(ni) += npos2;
	lengths_(ni) = npos2;    
      }    
    }
    else 
    {
      int Pos1 = pos1_;
      for ( int j = 0; j < ni; j++ )
      {
        if ( gaps_(j) < 0 ) Pos1 -= gaps_(j);
	Pos1 += lengths_(j);    
      }
      int npos1 = len1 - Pos1;
      Assert( npos1 >= 0 ); // XXX
      if ( npos1 == 0 ) --ni;
      else if ( npos1 >= gaps_(ni) )
      {
        lengths_(ni-1) += gaps_(ni);
	--ni;    
      }
      else
      {
        gaps_(ni) -= npos1;
	lengths_(ni) = npos1;    
      }    
    }    
  }
  
  nblocks_ = ni + 1;
  if ( gaps_(0) > 0 ) // shouldn't happen
  {
    pos2_ += gaps_(0);
    gaps_(0) = 0;    
  }
  else if ( gaps_(0) < 0 ) // shouldn't happen
  {
    pos1_ -= gaps_(0);
    gaps_(0) = 0;    
  }
  
 finish:
  
  // Make alignments proper if only off by one base at the beginning.
  // (This is really covering up a bug in some other piece of code.)
  
  if ( (pos1_ == 1 && pos2_ > 0) || (pos1_ > 0 && pos2_ == 1) )
  {
    --pos1_;
    --pos2_;
    ++lengths_(0);    
  }
  

  pos2_ += add_gap_plus;
  // Not sure about this one...
  //pos1_ += add_gap_minus;
  

  Bool zero_gap = False;
  for ( int i = 1; i < nblocks_; i++ )
    if ( gaps_(i) == 0 ) zero_gap = True;
  if (zero_gap) goto check_first_gap;    

}


void align::UnpackFrom( const packalign& p )
{
  p.Unpack( pos1_, pos2_, gaps_, lengths_, nblocks_ );    
}

// Proper tests an align to see if it is proper.  It does not flag "dead"
// alignments, those with no gaps or lengths.

Bool Proper( const align& a, int len1, int len2 )
{
  if ( a.Nblocks( ) == 0 ) return True;
  int pos1 = a.pos1( ), pos2 = a.pos2( ), Pos1 = a.Pos1( ), Pos2 = a.Pos2( );
  for ( int i = 1; i < a.Nblocks( ); i++ )
    if ( a.Gaps(i) == 0 ) return False;
  for ( int i = 0; i < a.Nblocks( ); i++ )
    if ( a.Lengths(i) < 0 || ( i < a.Nblocks( ) - 1 && a.Lengths(i) == 0 ) )
      return False;
  if ( a.Gaps(0) != 0 || Pos1 > len1 || Pos2 > len2 || pos1 < 0 || pos2 < 0 ||
       !( (pos1 == 0 && Pos2 == (int) len2)
          || (pos2 == 0 && Pos1 == (int) len1)
          || (pos1 == 0 && Pos1 == (int) len1)
          || (pos2 == 0 && Pos2 == (int) len2) ) )
  {
    return False;    
  }
  return True;    
}

// RequireProper tests an alignment_plus to see if it is proper.  It does not
// object to "dead" alignments, those with no gaps or lengths.

void RequireProper( const align& a, int id1, int id2, Bool rc2, 
                    const vecbasevector& EE, int test_no, Bool fatal )
{
  if ( Proper( a, EE[id1].size( ), EE[id2].size( ) ) ) return;
  int pos1 = a.pos1( ), pos2 = a.pos2( ), Pos1 = a.Pos1( ), Pos2 = a.Pos2( );
  PRINT2( id1, id2 );
  PRINT4( pos1, Pos1, pos2, Pos2 );
  PRINT2( EE[id1].size( ), EE[id2].size( ) );
  cout << "actual errors = " << flush;
  if ( !rc2 ) cout << ActualErrors( EE[id1], EE[id2], a ) << "\n";
  else 
  {
    basevector EErcid2 = EE[id2];
    EErcid2.ReverseComplement( );
    cout << ActualErrors( EE[id1], EErcid2, a ) << "\n";    
  }
  cout << "gaps/lengths:" << flush;
  for ( int i = 0; i < a.Nblocks( ); i++ )
    cout << " " << a.Gaps(i) << "/" << a.Lengths(i);
  cout << "\n";
  PRINT(test_no);
  if (fatal) ForceAssert( 0 == 1 );    
}

void RequireProper( const align& a, int id1, int id2, Bool rc2, 
                    const vec<int>& EE_length, int test_no, Bool fatal )
{
  if ( Proper( a, EE_length[id1], EE_length[id2] ) ) return;
  int pos1 = a.pos1( ), pos2 = a.pos2( ), Pos1 = a.Pos1( ), Pos2 = a.Pos2( );
  PRINT2( id1, id2 );
  PRINT4( pos1, Pos1, pos2, Pos2 );
  PRINT2( EE_length[id1], EE_length[id2] );
  cout << "gaps/lengths:" << flush;
  for ( int i = 0; i < a.Nblocks( ); i++ )
    cout << " " << a.Gaps(i) << "/" << a.Lengths(i);
  cout << "\n";
  PRINT(test_no);
  if (fatal) ForceAssert( 0 == 1 );    
}

int ActualErrors( const basevector& rd1, const basevector& rd2, 
                  const align& a, int mismatch_penalty, int gap_penalty )
{    int p1 = a.pos1( ), p2 = a.pos2( ), answer = 0;
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 )
          {    p2 += a.Gaps(j);
               answer += gap_penalty * a.Gaps(j);    }
          if ( a.Gaps(j) < 0 )
          {    p1 -= a.Gaps(j);
               answer -= gap_penalty * a.Gaps(j);    }
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    if ( rd1[p1] != rd2[p2] ) answer += mismatch_penalty;
               ++p1; ++p2;    }    }
     return answer;    }

int ActualErrors( const vec<char>& rd1, const basevector& rd2, 
                  const align& a, int mismatch_penalty, int gap_penalty )
{    int p1 = a.pos1( ), p2 = a.pos2( ), answer = 0;
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 )
          {    p2 += a.Gaps(j);
               answer += gap_penalty * a.Gaps(j);    }
          if ( a.Gaps(j) < 0 )
          {    p1 -= a.Gaps(j);
               answer -= gap_penalty * a.Gaps(j);    }
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    if ( toupper(rd1[p1]) != as_base(rd2[p2]) ) 
                    answer += mismatch_penalty;
               ++p1; ++p2;    }    }
     return answer;    }

int ActualErrors( const fastavector& rd1, const basevector& rd2, 
                  const align& a, int mismatch_penalty, int gap_penalty )
{    int p1 = a.pos1( ), p2 = a.pos2( ), answer = 0;
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 )
          {    p2 += a.Gaps(j);
               answer += gap_penalty * a.Gaps(j);    }
          if ( a.Gaps(j) < 0 )
          {    p1 -= a.Gaps(j);
               answer -= gap_penalty * a.Gaps(j);    }
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    if ( !GeneralizedBase::fromChar(rd1[p1])
                    .matches( GeneralizedBase::fromChar( as_base(rd2[p2]) ) ) )
               {    answer += mismatch_penalty;    }
               ++p1; ++p2;    }    }
     return answer;    }

int ActualErrors( const fastavector& rd1, const fastavector& rd2, 
                  const align& a, int mismatch_penalty, int gap_penalty )
{    int p1 = a.pos1( ), p2 = a.pos2( ), answer = 0;
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 )
          {    p2 += a.Gaps(j);
               answer += gap_penalty * a.Gaps(j);    }
          if ( a.Gaps(j) < 0 )
          {    p1 -= a.Gaps(j);
               answer -= gap_penalty * a.Gaps(j);    }
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    if ( !GeneralizedBase::fromChar(rd1[p1])
                    .matches( GeneralizedBase::fromChar(rd2[p2]) ) ) 
               {    answer += mismatch_penalty;    }
               ++p1; ++p2;    }    }
     return answer;    }

int ActualErrorsRc( const basevector& rd1, const basevector& rd2, 
                    const align& a, int mismatch_penalty, int gap_penalty )
{    int n2 = rd2.size( );
     int p1 = a.pos1( ), p2 = a.pos2( ), answer = 0;
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 )
          {    p2 += a.Gaps(j);
               answer += gap_penalty * a.Gaps(j);    }
          if ( a.Gaps(j) < 0 )
          {    p1 -= a.Gaps(j);
               answer -= gap_penalty * a.Gaps(j);    }
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    if ( rd1[p1] != 3 - rd2[n2-p2-1] ) answer += mismatch_penalty;
               ++p1; ++p2;    }    }
     return answer;    }

int ActualErrorsRc( const fastavector& rd1, const basevector& rd2, 
                    const align& a, int mismatch_penalty, int gap_penalty )
{    int n2 = rd2.size( );
     int p1 = a.pos1( ), p2 = a.pos2( ), answer = 0;
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 )
          {    p2 += a.Gaps(j);
               answer += gap_penalty * a.Gaps(j);    }
          if ( a.Gaps(j) < 0 )
          {    p1 -= a.Gaps(j);
               answer -= gap_penalty * a.Gaps(j);    }
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    if ( !GeneralizedBase::fromChar(rd1[p1]).complement( )
                    .matches( GeneralizedBase::fromChar( as_base(rd2[n2-p2-1]) ) ) ) 
               {    answer += mismatch_penalty;    }
               ++p1; ++p2;    }    }
     return answer;    }

int ActualErrorsRc( const fastavector& rd1, const fastavector& rd2, 
                    const align& a, int mismatch_penalty, int gap_penalty )
{    int n2 = rd2.size( );
     int p1 = a.pos1( ), p2 = a.pos2( ), answer = 0;
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 )
          {    p2 += a.Gaps(j);
               answer += gap_penalty * a.Gaps(j);    }
          if ( a.Gaps(j) < 0 )
          {    p1 -= a.Gaps(j);
               answer -= gap_penalty * a.Gaps(j);    }
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    if ( !GeneralizedBase::fromChar(rd1[p1]).complement( )
                    .matches( GeneralizedBase::fromChar(rd2[n2-p2-1]) ) ) 
               {    answer += mismatch_penalty;    }
               ++p1; ++p2;    }    }
     return answer;    }

template<class BASEVEC1, class BASEVEC2>
int ActualErrors( Bool rc, const BASEVEC1& rd1, const BASEVEC2& rd2, 
                  const align& a, int mismatch_penalty, int gap_penalty )
{
  if ( !rc ) return ActualErrors( rd1, rd2, a, mismatch_penalty, gap_penalty );
  else return ActualErrorsRc( rd1, rd2, a, mismatch_penalty, gap_penalty );    
}

template int ActualErrors( Bool rc, const basevector& rd1, const basevector& rd2, 
                  const align& a, int mismatch_penalty, int gap_penalty );
template int ActualErrors( Bool rc, const fastavector& rd1, const basevector& rd2, 
                  const align& a, int mismatch_penalty, int gap_penalty );
template int ActualErrors( Bool rc, const fastavector& rd1, const fastavector& rd2, 
                  const align& a, int mismatch_penalty, int gap_penalty );

int Bandwidth( align& a )
{
  const int add_to_bandwidth = 8; // heuristic
  int low = 0, high = 0, gap_total = 0;
  for ( int l = 0; l < a.Nblocks( ); l++ )
  {
    gap_total += a.Gaps(l);
    low = Min( low, gap_total );
    high = Max( high, gap_total );    
  }
  return Max( Abs(low), Abs(high) ) + add_to_bandwidth;    
}

void align::Read( istream& in, int& errors, int& id1, int& id2, Bool& rc )
{
  int pos1, pos2, gap, len, nblocks;
  BinRead( in, pos1 );
  Setpos1(pos1);
  BinRead( in, pos2 );
  Setpos2(pos2);
  BinRead( in, errors );
  BinRead( in, nblocks );
  SetNblocks(nblocks);
  for ( int i = 0; i < nblocks; i++ )
  {
    BinRead( in, gap );
    SetGap( i, gap );
    BinRead( in, len );
    SetLength( i, len );    
  }
  BinRead( in, id1 );
  BinRead( in, id2 );
  BinRead( in, rc );    
}

int CorrelatePositions( const align& a, int x1 )
{
  int pos1 = a.pos1( ), pos2 = a.pos2( );
  if ( x1 < pos1 ) return OffTheEnd;
  if ( x1 == pos1 ) return pos2;
  int nblocks = a.Nblocks( );
  const avector<int>& gaps = a.Gaps( );
  const avector<int>& lengths = a.Lengths( );
  for ( int j = 0; j < nblocks; j++ )
  {
    if ( gaps(j) > 0 ) pos2 += gaps(j);
    if ( gaps(j) < 0 )
      for ( int x = 0; x < -gaps(j); x++ )
	if ( x1 == pos1++ ) return AtGap;
    if ( x1 < pos1 + lengths(j) ) return pos2 + x1 - pos1;
    else
    {
      pos1 += lengths(j);
      pos2 += lengths(j);    
    }    
  }
  return OffTheEnd;    
}

void align::ReverseThis( int b1_len, int b2_len )
{
  vec<int> gaps_new, lengths_new;
  gaps_new.resize(0);
  lengths_new.resize(0);
  int n = Nblocks( );
  for ( int i = n-1; i >= 0; i-- )
  {
    if ( i == n-1 ) 
    {
      gaps_new.push_back(0);
      lengths_new.push_back( lengths_(n-1) );    
    }
    else 
    {
      gaps_new.push_back( gaps_(i+1) );
      lengths_new.push_back( lengths_(i) );    
    }    
  }
  if ( n > 0 && gaps_(0) != 0 ) 
  {
    gaps_new.push_back( gaps_(0) );
    lengths_new.push_back(0);    
  }
  int Pos1 = pos1_, Pos2 = pos2_;
  for ( int j = 0; j < n; j++ )
  {
    if ( gaps_(j) > 0 ) Pos2 += gaps_(j);
    else if ( gaps_(j) < 0 ) Pos1 -= gaps_(j);
    Pos1 += lengths_(j); 
    Pos2 += lengths_(j);    
  }
  pos1_ = b1_len - Pos1;
  pos2_ = b2_len - Pos2;
  SetNblocks( gaps_new.size( ) );
  for ( int i = 0; i < nblocks_; i++ )
  {
    gaps_(i) = gaps_new[i];
    lengths_(i) = lengths_new[i];    
  }    
}


void align::Flip( )
{
  swap( pos1_, pos2_ );
  for ( int i = 0; i < nblocks_; i++ )
    gaps_(i) = -gaps_(i);    
}


align align::TrimmedTo1(const int start, const int len) const {
  //make sure start and end are sensible!
  int endOn1 = start + len;
  int startOn1 = start;
  //if the trimming goes beyond the alignment, reduce it to alignment size.
  if (endOn1 >= Pos1()) {
    endOn1 = Pos1()-1; 
  }
  if (startOn1 < pos1()) {
    startOn1 = pos1(); 
  }
  align ret;
  ret.pos1_ = ret.pos2_ = ret.nblocks_ = 0;
  //if the trimming no longer makes sense, return an empty alignment.
  if (endOn1 <= startOn1) return ret;

  int p1 = pos1( ), p2 = pos2( );
  vec<int> gaps, blocks; //for the new align.
  bool started=false;
  bool done = false;
  int block=0, local=0;
  bool onGap=true;//we start on a gap of length 0.
  while (true) {
    //keep track of where we are
    int currentSize = onGap ? abs(Gaps(block)) : Lengths(block);

    if (local == currentSize) {
      //switching blocks: save the upcoming one and reset counters.
      if (started) {
	if (!onGap) {
	  gaps.push_back(Gaps(block+1));
	  //cout << "added gap " << gaps.size()-1 << "of size " << gaps.back() << endl;
        }
	else {
	  blocks.push_back(Lengths(block));
	  //cout << "added block " << blocks.size()-1 << "of size " << blocks.back() << endl;
        }
      }
      if (!onGap) {
        ++block; 
      }
      onGap = !onGap;
      local = 0;
      currentSize = onGap ? abs(Gaps(block)) : Lengths(block);
    }
    if (p1 >= startOn1 
	&& !started //start only once!
	&& !onGap) {
      //We do not want to start on a gap!
      //so we just wait until we are in a block to turn started to true.
      started = true;
      startOn1 = p1;
      gaps.push_back(0);
      //cout << "Added gap "<< gaps.size()-1<<" with size "<<gaps.back()<<endl;
      blocks.push_back(currentSize - local);
      //cout << "Added block " << blocks.size()-1<<" with size "<<blocks.back()<<endl;
      ret.pos1_ = p1;
      ret.pos2_ = p2;
      //cout << "starting at " << (onGap ? "gap " : "block ")
      //<< block << " position " << local << " p1 is " << p1 << endl
      //<< "currentSize " << currentSize << endl;
    } // end of if starting now

    //note that it is possible to come to the end and not have started if
    //the whole extent of b1 we are interested in is within an insertion.
    if (p1 == endOn1) {
      //note that it is possible to come to the end and not have started if
      //the whole extent of b1 we are interested in is within an insertion.
      if (started) {
	if (onGap) {
	  //remove the ending gap.
	  gaps.resize(gaps.size()-1);
        } else {
	  int bsize=(blocks.size()) == 1 ? endOn1 - startOn1 : local;
	  blocks.back() = bsize;
        }
      }
      //cout << "ending at " << (onGap ? "gap " : "block ")
      //<< block << " position " << local << " p1 is " << p1 << endl;
      break;
    }
    //cout << "at: " << (onGap ? "gap " : "block ")
    //<< block << " position " << local << " p1 is " << p1
    //<< " currentSize " << currentSize << endl;


    //advance local, p1 and p2 appropriately.
    if (onGap) {
      if (Gaps(block) > 0) {
	++p2;
      }
      else if (Gaps(block) < 0) {
	++p1;    
      }
    } else {
      ++p1;
      ++p2;
    }
    ++local;
  }
  //Deal with the situation where we are at the end of a gap:
  if (!blocks.empty() && blocks.back() == 0) {
    blocks.resize(blocks.size()-1);
    gaps.resize(gaps.size()-1);
  }
  if (started) {
    ret.pos1_ = ret.pos1_ - start;//move the start to the trimmed read.
    ret.nblocks_ = gaps.size();
    ret.gaps_ = avector<int>(gaps.begin(), gaps.end());
    ret.lengths_ = avector<int>(blocks.begin(), blocks.end());
  }
  //PRINT4(ret.pos1_, ret.pos2_, gaps, blocks);
  return ret;    
}  

void Trim1Together(const basevector & b1, const basevector & b2, 
		   const align & a, int startOn1, int len, 
		   basevector & trimmedb1, align & trimmeda) {
  trimmeda = a.TrimmedTo1(startOn1, len);
  trimmedb1.SetToSubOf(b1, startOn1, len);
}


int align::Errors( const basevector& rd1, const basevector& rd2 ) const {
  vec<int> errs = MutationsGap1Gap2( rd1, rd2 );
  return accumulate(errs.begin(), errs.end(), 0);
}

Bool align::Perfect( const basevector& rd1, const basevector& rd2 ) const 
{    if ( Nblocks( ) != 1 ) return False;
     if ( Gaps(0) != 0 ) return False; // would be weird
     int p1 = pos1( ), p2 = pos2( );
     for ( int x = 0; x < Lengths(0); x++ ) 
     {    if ( rd1[p1] != rd2[p2] ) return False;
          ++p1; ++p2;    }
     return True;    }

vector<int> align::MutationsGap1Gap2( const basevector& rd1, 
                                      const basevector& rd2 ) const {
  vector<int> answer(3, 0);
  int p1 = pos1( ), p2 = pos2( );
  for ( int j = 0; j < Nblocks( ); j++ ) {
    if ( Gaps(j) > 0 )  {
      answer[1] += Gaps(j);
      p2 += Gaps(j);    
    }
    if ( Gaps(j) < 0 ) {
      answer[2] -= Gaps(j);
      p1 -= Gaps(j);    
    }
    for ( int x = 0; x < Lengths(j); x++ ) {
      if ( rd1[p1] != rd2[p2] ) ++answer[0];
      ++p1; 
      ++p2;    
    }    
  }
  return answer;    
}

int align::PosOn1(int on2) const {
  if (on2 < pos2() || on2 > Pos2()) return -1;
  int p1 = pos1(), p2 = pos2();
  for ( int j = 0; j < Nblocks( ); j++ ) {
    if ( Gaps(j) > 0 )  {
      p2 += Gaps(j);
    }
    if ( Gaps(j) < 0 ) {
      p1 -= Gaps(j);
    }
    for ( int x = 0; x < Lengths(j); x++ ) {
      if ( p2 == on2 ) return p1;
      if (p2 > on2) return -1;
      ++p1;
      ++p2;
    }
  }
  return -1;
}

int align::PosOn2(int on1) const {
  if (on1 < pos1() || on1 > Pos1()) return -1;
  int p1 = pos1(), p2 = pos2();
  for ( int j = 0; j < Nblocks( ); j++ ) {
    if ( Gaps(j) > 0 )  {
      p2 += Gaps(j);
    }
    if ( Gaps(j) < 0 ) {
      p1 -= Gaps(j);
    }
    for ( int x = 0; x < Lengths(j); x++ ) {
      if ( p1 == on1 ) return p2;
      if (p1 > on1) return -1;
      ++p1;
      ++p2;
    }
  }
  return -1;
}

pair<int, int> align::Gap1Gap2( ) const {
  pair<int, int> ret(0,0);
  for ( int j = 0; j < Nblocks( ); j++ ) {
    if ( Gaps(j) > 0 )  {
      ret.first += Gaps(j);  
    }
    else if ( Gaps(j) < 0 ) {
      ret.second -= Gaps(j);   
    }
  }
  return ret;    
}



void
align::Sync_to_TACG( const basevector & seq1,
                     const basevector & seq2,
                     Bool  isRC )
{
  // Adjust the align object to synchronize with the 454-cycles of TACG.
  //
  // Upon return from this method, the members  nblocks_, gaps_  and  lengths_
  // will be adjusted so that the extents of the alignment on both sequences
  // are the same as before, but there is no mismatch within each block;
  // instead, mismatches will manifest themselves as "extra" indels
  // in the align object, in sync with the flow order TACG on  seq2.
  //
  // If isRC is False, the start positions (pos1(), pos2()) of the two sequences
  // are fixed but the end positions (Pos1(), Pos2()) may change.
  // If isRC is True,  the end positions (Pos1(), Pos2()) of the two sequences
  // are fixed but the start positions (pos1(), pos2()) may change.
  //

  const unsigned char flow_order[4] = { 3, 0, 1, 2 };
  // Via the dictionary  0 = A, 1 = C, 2 = G, 3 = T
  // this corresponds to the flow order TACG.

  const int start1 = pos1();
  const int start2 = pos2();
  const int end1 = Pos1();
  const int end2 = Pos2();

  const int max_flows = 4 * Max(end1 - start1, end2 - start2);

  vec<int> flow1(max_flows,0);
  vec<int> flow2(max_flows,0);
  int flow_index = 0;

  // convert the two subsequences in alignment into flow multiplicities
  if ( isRC == False )
  {
    // If  isRC == False, then the synchronization of cycles starts with
    // the beginning of  seq1  and the beginning of  seq2;

    int base_index1 = start1;
    int base_index2 = start2;

    while ( base_index1 < end1 ||
            base_index2 < end2 )
    {
      const unsigned char flow_base = flow_order[ flow_index % 4 ];

      unsigned int freq = 0;
      while ( base_index1 < end1 &&
              seq1[base_index1] == flow_base )
      {
	base_index1++;
	freq++;
      }
      flow1[ flow_index ] = freq;

      freq = 0;
      while ( base_index2 < end2 &&
              seq2[base_index2] == flow_base )
      {
	base_index2++;
	freq++;
      }
      flow2[ flow_index ] = freq;

      flow_index++;
    }
  }      
  else // if ( isRC == True )
  {
    // If  isRC == True,  then the synchronization of cycles starts with
    // the end       of  seq1  and the end       of  seq2;

    int base_index1 = end1 - 1;
    int base_index2 = end2 - 1;

    while ( base_index1 >= start1 ||
            base_index2 >= start2 )
    {
      const unsigned char flow_base = flow_order[ flow_index % 4 ];

      unsigned int freq = 0;
      while ( base_index1 >= start1 &&
              seq1[base_index1] == flow_base )
      {
	base_index1--;
	freq++;
      }
      flow1[ flow_index ] = freq;

      freq = 0;
      while ( base_index2 >= start2 &&
              seq2[base_index2] == flow_base )
      {
	base_index2--;
	freq++;
      }
      flow2[ flow_index ] = freq;

      flow_index++;
    }
  }

  const int num_flows = flow_index;
  flow1.resize(num_flows);
  flow2.resize(num_flows);


  // compute the expanded vectors of length/gap info;
  vec<int> gaps(num_flows+1, 0);
  vec<int> lengths(num_flows+1, 0);

  for (flow_index = 0; flow_index < num_flows; flow_index++)
  {
    lengths[flow_index] = Min( flow1[flow_index],
                               flow2[flow_index] );

    if ( flow1[flow_index] > lengths[flow_index] )
    {
      gaps[flow_index] = - (flow1[flow_index] - lengths[flow_index]);
    }
    else if (flow2[flow_index] > lengths[flow_index] )
    {
      gaps[flow_index] = + (flow2[flow_index] - lengths[flow_index]);
    }
  }
  // lengths[i]  now contains the number (>= 0) of common bases
  //             in the i-th flow;
  // gaps[i]     now contains the number of bases in the i-th flow
  //             which flow2 has in excess of flow1


  // compactify length/gap info
  int block_index = 0;

  flow_index = 0;
  while (flow_index < num_flows)
  {
    // advance through a contiguous length-block (i.e. no gaps)
    int curr_length = lengths[flow_index];
    while ( flow_index < num_flows &&
            gaps[flow_index] == 0 )
    {
      flow_index++;
      curr_length += lengths[flow_index];
    }

    int curr_gap = gaps[flow_index];
    flow_index++;

    // advance through a contiguous gap-block
    // (i.e. no lengths, and gaps are of the same sign)
    while ( flow_index < num_flows &&
            lengths[flow_index] == 0 &&
            gaps[flow_index] * curr_gap >= 0 )
    {
      curr_gap += gaps[flow_index];
      flow_index++;
    }

    // curr_length  and  curr_gap  now contain the length and gap sizes
    // for the current block (indexed by block_index);
    // we overwrite the earlier part of the lengths/gaps vector
    // to compactify

    lengths[block_index] = curr_length;
    gaps[block_index] = curr_gap;

    // advance to a new block
    block_index++;
  }

  // lengths[i]  and  gaps[i]  for i = 0,...,block_index-1
  // now contain  the compactified length/gap info;
  // we now copy that over to the private members of the align class

  const int num_blocks = block_index;
  SetNblocks(num_blocks);
  SetGap(0,0);

  if ( isRC == False )
  {
    SetLength(0, lengths[0]);
    for (int i = 1; i < num_blocks; i++)
    {
      SetGap(i, gaps[i-1]);
      SetLength(i, lengths[i]);
    }
  }
  else // if ( isRC == True )
  {
    // also compute the new starting positions
    int newstart1 = end1;
    int newstart2 = end2;

    for (int i = 0; i < num_blocks-1; i++)
    {
      SetLength(num_blocks-1-i, lengths[i]);
      SetGap(num_blocks-1-i, gaps[i]);

      newstart1 -= lengths[i];
      newstart2 -= lengths[i];

      if (gaps[i] > 0)
      {
	newstart2 -= gaps[i];
      }
      else if (gaps[i] < 0)
      {

	newstart1 += gaps[i];
      
      }
    }
    SetLength(0, lengths[num_blocks-1]);
    newstart1 -= lengths[num_blocks-1];
    newstart2 -= lengths[num_blocks-1];

    Setpos1(newstart1);
    Setpos2(newstart2);
  }
  // done!
}



int align::Mutations( const basevector& rd1, const basevector& rd2,
                      const qualvector& q1, int min_score ) const
{
  int answer = 0, j, p1 = pos1( ), p2 = pos2( );
  for ( j = 0; j < Nblocks( ); j++ )
  {
    if ( Gaps(j) > 0 ) p2 += Gaps(j);
    if ( Gaps(j) < 0 ) p1 -= Gaps(j);
    for ( int x = 0; x < Lengths(j); x++ )
    {
      if ( rd1[p1] != rd2[p2] && q1[p1] >= min_score ) ++answer;
      ++p1; 
      ++p2;    
    }    
  }
  return answer;    
}

void align::PrintMutations( const basevector& rd1, const basevector& rd2, ostream& log) const
{
  int answer = 0, j, p1 = pos1( ), p2 = pos2( );
  for ( j = 0; j < Nblocks( ); j++ )
  {
    if ( Gaps(j) > 0 ) {
      log << p1 << " insertion " << p2 << " ";
      for (int i = 0; i < Gaps(j); ++i)
	log << Base::val2Char(rd2[p2+i]);
      log << endl;
      p2 += Gaps(j);
    }
    if ( Gaps(j) < 0 ) {
      log << p1 << " ";
      for (int i = 0; i < -Gaps(j); ++i)
	log << Base::val2Char(rd1[p1+i]);
      log << " "  << p2 << " deletion" << endl;
      p1 -= Gaps(j);
    }
    for ( int x = 0; x < Lengths(j); x++ )
    {
      if ( rd1[p1] != rd2[p2] ) {
	// base mismatch
	log << p1 << " " << Base::val2Char(rd1[p1])
	    << " " << p2 << " " << Base::val2Char(rd2[p2]) << endl;
      }
      ++p1; 
      ++p2;    
    }    
  }
}

int align::MatchingBases( const basevector& rd1, const basevector& rd2 )
{
  int answer = 0, j, p1 = pos1( ), p2 = pos2( );
  for ( j = 0; j < Nblocks( ); j++ )
  {
    if ( Gaps(j) > 0 ) p2 += Gaps(j);
    if ( Gaps(j) < 0 ) p1 -= Gaps(j);
    for ( int x = 0; x < Lengths(j); x++ )
    {
      if ( rd1[p1] == rd2[p2] ) ++answer;
      ++p1; 
      ++p2;    
    }    
  }
  return answer;    
}

void align::PerfectIntervals1( const basevector& rd1, const basevector& rd2,
     vec<ho_interval>& perfs ) const
{    perfs.clear( );
     int p1 = pos1( ), p2 = pos2( );
     for ( int j = 0; j < Nblocks( ); j++ )
     {    if ( Gaps(j) > 0 ) p2 += Gaps(j);
          if ( Gaps(j) < 0 ) p1 -= Gaps(j);
          int last = p1 - 1;
          for ( int x = 0; x < Lengths(j); x++ )
          {    if ( rd1[p1] != rd2[p2] )
               {    if ( p1 - (last+1) > 0 )
                         perfs.push_back( ho_interval( last+1, p1 ) );
	            last = p1;    }
               ++p1; 
               ++p2;    }
          if ( p1 - (last+1) > 0 ) 
               perfs.push_back( ho_interval( last+1, p1 ) );    }    }

void align::PerfectIntervals2( const basevector& rd1, const basevector& rd2,
     vec<ho_interval>& perfs ) const
{    perfs.clear( );
     int p1 = pos1( ), p2 = pos2( );
     for ( int j = 0; j < Nblocks( ); j++ )
     {    if ( Gaps(j) > 0 ) p2 += Gaps(j);
          if ( Gaps(j) < 0 ) p1 -= Gaps(j);
          int last = p2 - 1;
          for ( int x = 0; x < Lengths(j); x++ )
          {    if ( rd1[p1] != rd2[p2] )
               {    if ( p2 - (last+1) > 0 )
                         perfs.push_back( ho_interval( last+1, p2 ) );
	            last = p2;    }
               ++p1; ++p2;    }
          if ( p2 - (last+1) > 0 ) 
               perfs.push_back( ho_interval( last+1, p2 ) );    }    }

void align::PerfectIntervals2( const fastavector& rd1, const fastavector& rd2,
     vec<ho_interval>& perfs ) const
{    perfs.clear( );
     int p1 = pos1( ), p2 = pos2( );
     for ( int j = 0; j < Nblocks( ); j++ )
     {    if ( Gaps(j) > 0 ) p2 += Gaps(j);
          if ( Gaps(j) < 0 ) p1 -= Gaps(j);
          int last = p2 - 1;
          for ( int x = 0; x < Lengths(j); x++ )
          {    if ( !GeneralizedBase::fromChar(rd1[p1])
                    .matches( GeneralizedBase::fromChar(rd2[p2]) ) ) 
               {    if ( p2 - (last+1) > 0 )
                         perfs.push_back( ho_interval( last+1, p2 ) );
	            last = p2;    }
               ++p1; ++p2;    }
          if ( p2 - (last+1) > 0 ) 
               perfs.push_back( ho_interval( last+1, p2 ) );    }    }

int align::Indels( const basevector& rd1, const basevector& rd2,
                   const qualvector& q1, int min_score ) const
{
  int answer = 0, p1 = pos1( ), p2 = pos2( );
  for ( int j = 0; j < Nblocks( ); j++ )
  {
    if ( Gaps(j) > 0 ) 
    {
      if ( p1+1 < (int) q1.size( ) && q1[p1] >= min_score 
           && q1[p1+1] >= min_score ) answer += Gaps(j);
      p2 += Gaps(j);    
    }
    if ( Gaps(j) < 0 ) 
    {
      if ( -Gaps(j) + p1 < (int) q1.size( ) )
      {
        int k;
        for ( k = 0; k <= -Gaps(j); k++ )
          if ( q1[ p1 + k ] < min_score ) break;
        if ( k == -Gaps(j) + 1 ) answer -= Gaps(j);    
      }
      p1 -= Gaps(j);    
    }
    p1 += Lengths(j);
    p2 += Lengths(j);    
  }
  return answer;    
}

vector<int> align::MutationsGap1Gap2( const basevector& rd1, 
                                      int from1, int to1, 
                                      const basevector& rd2, 
                                      int from2, int to2 ) const
{
  vector<int> answer(3, 0);
  int p1 = pos1( ), p2 = pos2( );
  for ( int j = 0; j < Nblocks( ); ++j ) 
  {
    if ( Gaps(j) > 0 ) 
    {
      if ( p1 >= from1 && p1 < to1 && p2 >= from2 && p2 < to2 )
        answer[1] += Gaps(j);
      p2 += Gaps(j);    
    }
    if ( Gaps(j) < 0 ) 
    {
      if ( p1 >= from1 && p1 < to1 && p2 >= from2 && p2 < to2 )
	answer[2] -= Gaps(j);
      p1 -= Gaps(j);    
    }
    for ( int x = 0; x < Lengths(j); ++x ) 
    {
      if ( p1 >= from1 && p1 < to1 && p2 >= from2 && p2 < to2 &&
           rd1[p1] != rd2[p2] )
        ++answer[0];
      ++p1;
      ++p2;    
    }    
  }
  return answer;    
}

void align::Write( ostream& out, int id1, int id2, Bool rc, int errors )
{
  BinWrite( out, pos1_ );
  BinWrite( out, pos2_ );
  BinWrite( out, errors );
  int nblocks = Nblocks( );
  BinWrite( out, nblocks );
  int gap, length;
  for ( int i = 0; i < nblocks; i++ )
  {
    gap = Gaps(i);
    length = Lengths(i);
    BinWrite( out, gap );
    BinWrite( out, length );    
  }
  BinWrite( out, id1 );
  BinWrite( out, id2 );
  BinWrite( out, rc );    
}

ostream & operator<<(ostream & os, const align & a) {
  os << "startOn1: " << a.pos1() << ", startOn2: " << a.pos2()
     << ", Nblocks: " << a.Nblocks() << endl;
  return os;
}

BINARY3_DEF(placement_mark);
template class SmallVec< placement_mark, MempoolAllocator<placement_mark> >;
template class OuterVec<PlacementMarkVec>;
