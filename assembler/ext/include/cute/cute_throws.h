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
 * Copyright 2007-2009 Peter Sommerlad
 *
 *********************************************************************************/

#ifndef CUTE_THROWS_H_
#define CUTE_THROWS_H_
#include "cute_base.h"

// should we allow arbitrary code and remove the parentheses around the macro expansion?
// not now, strange compilation side-effects might result.
namespace cute {
	namespace do_not_use_this_namespace {
		struct assert_throws_failure_exception {
			struct cute::test_failure original;
			assert_throws_failure_exception(std::string const &r,char const *f, int line):
					original(r,f, line){}
		};
	}
}
#define ASSERT_THROWSM(msg,code,exc) \
	do { \
		try { \
			{ code ; } \
			throw cute::do_not_use_this_namespace::assert_throws_failure_exception(#msg,__FILE__,__LINE__); \
		} catch(exc &){ \
		} catch(cute::do_not_use_this_namespace::assert_throws_failure_exception &atf){throw atf.original;} \
	} while(0)
#define ASSERT_THROWS(code,exc) ASSERT_THROWSM(" expecting " #code " to throw " #exc,code,exc)

#endif /*CUTE_THROWS_H_*/
