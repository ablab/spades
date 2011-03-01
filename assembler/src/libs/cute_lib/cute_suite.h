/*********************************************************************************
 * This file is part of CUTE.
 *
 * CUTE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CUTE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with CUTE.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2006 Peter Sommerlad
 *
 *********************************************************************************/
#ifndef CUTE_SUITE_H_
#define CUTE_SUITE_H_
#include "cute_test.h"
#include <vector>
namespace cute {
	typedef std::vector<test> suite;
	// convenience operator for appending to suites, might not be right
	// can use boost/assign.hpp instead...
	inline
	suite &operator+=(suite &left, suite const &right){
		left.insert(left.end(),right.begin(),right.end());
		return left;
	}
	inline
	suite &operator+=(suite &left, test const &right){
		left.push_back(right);
		return left;
	}
}
#endif /*CUTE_SUITE_H_*/
