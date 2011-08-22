/*
 * qualTest.hpp
 *
 *  Created on: 03.03.2011
 *      Author: vyahhi
 */

#ifndef QUALTEST_HPP_
#define QUALTEST_HPP_

#include "quality.hpp"

void TestQuality() {
	Quality q("0123456789");
	ASSERT_EQUAL('0', q[0]);
	ASSERT_EQUAL('6', q[6]);
	ASSERT_EQUAL('9', q[9]);
}

cute::suite QualitySuite(){
	cute::suite s;
	s.push_back(CUTE(TestQuality));
	return s;
}

#endif /* QUALTEST_HPP_ */
