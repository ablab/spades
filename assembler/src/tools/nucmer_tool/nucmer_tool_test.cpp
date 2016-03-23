//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * nucmer_tool_test.cpp
 *
 *  Created on: 22.08.2011
 *
 */
#include "nucmer_coords_tool.hpp"
#include <iostream>
#include <algorithm>
using namespace std;


int main() {
    std::vector<NucmerAllign> test_vect = ReadCoordsFile("test.coords");
    std::cout<<test_vect.size()<<endl;
    std::sort(test_vect.begin(), test_vect.end(), MyLess);
    NucmerAllign last = *test_vect.begin();
    int j = 0;    
    for (int i = 1; i<4700; i++){
        bool intersect = false;
        while((j<test_vect.size())&&(test_vect[j].inReferenceEnd < (i - 1)*1000)){    
//            printf("skip %i %i\n", test_vect[j].inReferenceStart, test_vect[j].inReferenceEnd);
            j++;
        }
//        printf("check %i %i\n", test_vect[j].inReferenceStart, test_vect[j].inReferenceEnd);
        if ( (test_vect[j].inReferenceStart < (i)*1000) && ((test_vect[j].inReferenceEnd > (i - 1)*1000)) )
            printf("%i 1\n", i);
        else
            printf("%i 0\n" , i);

    }
//    for (auto i = test_vect.begin(); i<test_vect.end(); ++i){
//        if (!Redundant(last, *i) || Redundant(*i, last)){
//            printf("%8i %8i | %8i %8i | %8i %8i | %8.2f | %s", i->inReferenceStart, i->inReferenceEnd
//                                                     , i->inContigStart, i->inContigEnd
//                                                     , i->inReferenceLength, i->inContigLength
//                                                     , i->ReferenceQuality, i->contigName.c_str());
//            last = *i;
//        }
//    }
    return 0;
}
