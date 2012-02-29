///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "system/Types.h"
#include "system/StaticAssert.h"

#include <iostream>

using std::cerr;

void TestTypes( )
{
  STATIC_ASSERT_M(sizeof(char) == 1, bad_char_size);
  STATIC_ASSERT_M(sizeof(short) == 2, bad_short_size);
  STATIC_ASSERT_M(sizeof(int) == 4, bad_int_size);
  STATIC_ASSERT_M(sizeof(longlong) == 8, bad_longlong_size);
}
