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

#ifndef CUTE_TESTMEMBER_H_
#define CUTE_TESTMEMBER_H_
#include "cute_test.h"
#if ! defined(USE_TR1) && ! defined(USE_STD0X)
#include <boost/bind.hpp>
#endif
namespace cute {
	template <typename TestClass>
	test makeMemberFunctionTest(TestClass &t,void (TestClass::*fun)(),char const *name){
		return test(boost_or_tr1::bind(fun,boost_or_tr1::ref(t)),demangle(typeid(TestClass).name())+"::"+name);
	}
	template <typename TestClass>
	test makeMemberFunctionTest(TestClass const &t,void (TestClass::*fun)()const,char const *name){
		return test(boost_or_tr1::bind(fun,boost_or_tr1::cref(t)),demangle(typeid(TestClass).name())+"::"+name);
	}
	template <typename TestClass,typename MemFun>
	struct incarnate_for_member_function {
		MemFun memfun;
		incarnate_for_member_function(MemFun f):memfun(f){}
		void operator()(){
			TestClass t;
			(t.*memfun)();
		}
	};
	template <typename TestClass, typename MemFun>
	test makeSimpleMemberFunctionTest(MemFun fun,char const *name){
		return test(incarnate_for_member_function<TestClass,MemFun>(fun),demangle(typeid(TestClass).name())+"::"+name);
	}
	template <typename TestClass,typename MemFun, typename Context>
	struct incarnate_for_member_function_with_context_object {
		MemFun memfun;
		Context context;
		incarnate_for_member_function_with_context_object(MemFun f,Context c)
		:memfun(f),context(c){}
		void operator()(){
			TestClass t(context);
			(t.*memfun)();
		}
	};
	template <typename TestClass, typename MemFun, typename Context>
	test makeMemberFunctionTestWithContext(Context c,MemFun fun,char const *name){
		return test(incarnate_for_member_function_with_context_object<TestClass,MemFun,Context>(fun,c),demangle(typeid(TestClass).name())+"::"+name);
	}
}
#define CUTE_MEMFUN(testobject,TestClass,MemberFunctionName) \
	cute::makeMemberFunctionTest(testobject,\
		&TestClass::MemberFunctionName,\
		#MemberFunctionName)
#define CUTE_SMEMFUN(TestClass,MemberFunctionName) \
	cute::makeSimpleMemberFunctionTest<TestClass>(\
		&TestClass::MemberFunctionName,\
		#MemberFunctionName)
#define CUTE_CONTEXT_MEMFUN(context_object,TestClass,MemberFunctionName) \
	cute::makeMemberFunctionTestWithContext<TestClass>(\
		context_object,\
		&TestClass::MemberFunctionName,\
		#MemberFunctionName)

#endif /*CUTE_TESTMEMBER_H_*/
