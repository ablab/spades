/*
 * nuclTest.hpp
 *
 *  Created on: 02.03.2011
 *      Author: vyahhi
 */

#include "cute.h"
#include "nucl.hpp"

void TestNucl() {
	ASSERT_EQUAL('A', nucl(0));
	ASSERT_EQUAL('C', nucl(1));
	ASSERT_EQUAL('G', nucl(2));
	ASSERT_EQUAL('T', nucl(3));
	ASSERT_EQUAL(0, dignucl('A'));
	ASSERT_EQUAL(1, dignucl('C'));
	ASSERT_EQUAL(2, dignucl('G'));
	ASSERT_EQUAL(3, dignucl('T'));
	ASSERT_EQUAL(3, complement(0));
	ASSERT_EQUAL(2, complement(1));
	ASSERT_EQUAL(1, complement(2));
	ASSERT_EQUAL(0, complement(3));
	ASSERT(is_nucl('A'));
	ASSERT(is_nucl('C'));
	ASSERT(is_nucl('G'));
	ASSERT(is_nucl('T'));
	ASSERT(!is_nucl(0));
	ASSERT(!is_nucl(1));
	ASSERT(!is_nucl(2));
	ASSERT(!is_nucl(3));
	ASSERT(!is_nucl('0'));
	ASSERT(!is_nucl('1'));
}

cute::suite NuclSuite(){
	cute::suite s;
	s.push_back(CUTE(TestNucl));
	return s;
}
