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
#ifndef CUTE_TEST_INCARNATE_H_
#define CUTE_TEST_INCARNATE_H_
#include "cute_test.h"
// idea blatantly stolen from Aeryn
namespace cute {
	template <typename TestFunctor>
	struct test_incarnate {
		void operator()(){
			TestFunctor()();
		}
	};
	// TODO: check if there are problems with references.
	template <typename TestFunctor, typename ContextObject>
	struct test_incarnate_with_context {
		test_incarnate_with_context(ContextObject context):theContext(context)
		{}
		void operator()(){
			TestFunctor t(theContext);// wouldn't create temporary to call with ()()
			t();
		}
		ContextObject theContext;
	};
	template <typename TestFunctor,typename ContextObject>
	test make_incarnate_with_context(ContextObject obj){
		return test(test_incarnate_with_context<TestFunctor,ContextObject>(obj),demangle(typeid(TestFunctor).name()));
	}
}
#define CUTE_INCARNATE(TestFunctor) cute::test(cute::test_incarnate<TestFunctor>(),cute::demangle(typeid(TestFunctor).name()))
#define CUTE_INCARNATE_WITH_CONTEXT(TestFunctor,contextObject) cute::make_incarnate_with_context<TestFunctor>(contextObject)
#endif /*CUTE_TEST_INCARNATE_H_*/
