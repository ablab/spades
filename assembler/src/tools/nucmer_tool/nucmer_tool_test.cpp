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
	for (auto i = test_vect.begin(); i<test_vect.end(); ++i){
		if (!Redundant(last, *i) || Redundant(*i, last)){
			printf("%8i %8i | %8i %8i | %8i %8i | %8.2f | %s", i->inReferenceStart, i->inReferenceEnd
													 , i->inContigStart, i->inContigEnd
													 , i->inReferenceLength, i->inContigLength
													 , i->ReferenceQuality, i->contigName.c_str());
			last = *i;
		}
	}
	return 0;
}
