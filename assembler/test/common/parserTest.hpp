/*
 * parserTest.hpp
 *
 *  Created on: 02.03.2011
 *      Author: vyahhi
 */

#include "cute.h"
#include "parser.hpp"

void TestParserEmpty() {
	//
}

cute::suite ParserSuite(){
	cute::suite s;
	s.push_back(CUTE(TestParserEmpty));
	return s;
}
