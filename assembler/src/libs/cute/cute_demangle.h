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
 * Copyright 2009 Peter Sommerlad
 *
 *********************************************************************************/

#ifndef CUTE_DEMANGLE_H_
#define CUTE_DEMANGLE_H_
#include <string>
// needs adaptation for different compilers
// dependency to demangle is a given,
// otherwise we have to use macros everywhere
#ifdef __GNUG__
#include <cxxabi.h> // __cxa_demangle
#include <cstdlib> // ::free() 
namespace cute {

inline std::string demangle(char const *name){
	if (!name) return "unknown";
	char *toBeFreed = abi::__cxa_demangle(name,0,0,0);
	std::string result(toBeFreed?toBeFreed:name);
	::free(toBeFreed);
	return result;
}
}
#else
namespace cute {
// this default works reasonably with MSVC71 and 8, hopefully for others as well
inline std::string demangle(char const *name){
	return std::string(name?name:"unknown");
}
}
#endif

#endif /* CUTE_DEMANGLE_H_ */
