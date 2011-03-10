/*
 * nuclTest.hpp
 *
 *  Created on: 02.03.2011
 *      Author: vyahhi
 */

#include "cute.h"
#include "nucl.hpp"

void TestIt() {
	ASSERT_EQUAL('A', nucl(0));
	ASSERT_EQUAL('C', nucl(1));
	ASSERT_EQUAL('G', nucl(2));
	ASSERT_EQUAL('T', nucl(3));
	ASSERT_EQUAL(0, unnucl('A'));
	ASSERT_EQUAL(1, unnucl('C'));
	ASSERT_EQUAL(2, unnucl('G'));
	ASSERT_EQUAL(3, unnucl('T'));
	ASSERT_EQUAL(3, complement(0));
	ASSERT_EQUAL(2, complement(1));
	ASSERT_EQUAL(1, complement(2));
	ASSERT_EQUAL(0, complement(3));
}

cute::suite NuclSuite(){
	cute::suite s;
	s.push_back(CUTE(TestIt));
	return s;
}
