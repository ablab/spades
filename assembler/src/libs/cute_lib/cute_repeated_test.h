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

#ifndef CUTE_REPEATED_TEST_H_
#define CUTE_REPEATED_TEST_H_
#include "cute_test.h"
namespace cute{
	struct repeated_test {
		repeated_test(test const &t,unsigned int n):theTest(t),repetitions(n){}
		void operator()(){
			for (unsigned int i=0;i<repetitions;++i){
				theTest();
			}
		}
		test theTest;
		const unsigned int repetitions;
	};
}
#define CUTE_REPEAT(aTestFunction,n) cute::test(cute::repeated_test((aTestFunction),(n)),#aTestFunction " " #n " times repeated")
#define CUTE_REPEAT_TEST(aTestObject,n) cute::test(cute::repeated_test((aTestObject),(n)),aTestObject.name()+" " #n " times repeated")
#endif /*CUTE_REPEATED_TEST_H_*/
