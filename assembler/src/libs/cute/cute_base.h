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

#ifndef CUTE_BASE_H_
#define CUTE_BASE_H_
#include <string>
namespace cute{
	struct test_failure {
		std::string reason;
		std::string filename;
		int lineno;

		test_failure(std::string const &r,char const *f, int line)
		:reason(r),filename(f),lineno(line)
		{ 	}
		char const * what() const { return reason.c_str(); }
	};
}
#define ASSERTM(msg,cond) if (!(cond)) throw cute::test_failure((msg),__FILE__,__LINE__)
#define ASSERT(cond) ASSERTM(#cond,cond)
#define FAIL() ASSERTM("FAIL()",false)
#define FAILM(msg) ASSERTM(msg,false)
#endif /*CUTE_BASE_H_*/
