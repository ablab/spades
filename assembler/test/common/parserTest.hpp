/*
 * parserTest.hpp
 *
 *  Created on: 02.03.2011
 *      Author: vyahhi
 */

#include "cute.h"
#include "parser.hpp"

void TestEmpty() {
	//
}

cute::suite ParserSuite(){
	cute::suite s;
	s.push_back(CUTE(TestEmpty));
	return s;
}
