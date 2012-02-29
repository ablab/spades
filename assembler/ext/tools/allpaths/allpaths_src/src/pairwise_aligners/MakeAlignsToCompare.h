// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef MAKEALIGNSTOCOMPARE_H
#define MAKEALIGNSTOCOMPARE_H

#include "system/Assert.h"

// ================================================================================
//
// In MakeAligns, the argument which_to_compare should be:
//
//      to_compare(ALL_VS_ALL),
//      to_compare(FIRST_VS_SECOND, N0), or
//      to_compare(FIRST_VS_ALL, N0),
//
// where N0 is a nonnegative integer.
//
// In the ALL_VS_ALL case, all sequences are compared to each other.
//
// In the FIRST_VS_SECOND case, all sequences EE[i] with i < N0 are compared to all
// sequences EE[i] with i >= N0.
//
// In the FIRST_VS_ALL case, all sequences EE[i] with i < N0 are compared to all
// other sequences.
//
// ================================================================================

enum to_compare_type { ALL_VS_ALL, FIRST_VS_SECOND, FIRST_VS_ALL };

class to_compare {

     public:

     to_compare_type compare_type;
     int N0;

     to_compare( to_compare_type compare_type_arg,
		 int N0_arg = 0 )
       : compare_type(compare_type_arg),
	 N0(N0_arg)
     {
       ForceAssert( compare_type != ALL_VS_ALL || N0 == 0 );   
     }

};

#endif
