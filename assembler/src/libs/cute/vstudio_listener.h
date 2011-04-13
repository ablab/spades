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
 * Copyright 2006-2009 Peter Sommerlad
 *
 *********************************************************************************/

#ifndef VSTUDIO_LISTENER_H
#define VSTUDIO_LISTENER_H
// Windows listener for debug mode: allows selection of assert failing source line
// TODO: vstudio_listener is broken for VS later than 2003, because OutputDebugString no longer works as before
#ifndef __GNUG__
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <sstream>
#include <iostream>
namespace cute{
	class vstudio_listener
	{
	public:
		void begin(suite const &t,char const *info){
		}
		void end(suite const &t, char const *info){
		}
		void start(test const &t){
		}
		void success(test const &t, char const *msg){
			std::cerr <<  t.name() <<" " << msg<< std::endl;
		}
		void failure(test const &t,test_failure const &e){
			std::ostringstream out;
			out << e.filename << "(" << e.lineno << ") : testcase failed: " <<e.reason << " in " << t.name()<< std::endl;
			OutputDebugString(out.str().c_str());
			std::cerr << out.str() << std::flush;
		}
		void error(test const &t, char const *what){
			std::ostringstream out;
			out << what << " in " << t.name() << std::endl;
			OutputDebugString(out.str().c_str());
			std::cerr << out.str() << std::flush;
		}
	};
}
#else
// cheat for gnu use ostream_listener instead, so the type is defined
#include "ostream_listener.h"
namespace cute{
	typedef cute::ostream_listener vstudio_listener;
}
#endif
#endif
