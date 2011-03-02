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
 * Copyright 2007 Peter Sommerlad
 *
 *********************************************************************************/

#ifndef CUTE_COUNTING_LISTENER_H_
#define CUTE_COUNTING_LISTENER_H_
#include "cute_listener.h"
namespace cute{
	template <typename Listener=null_listener>
	struct counting_listener:Listener{
		counting_listener()
		:Listener()
		,numberOfTests(0),successfulTests(0),failedTests(0),errors(0),numberOfSuites(0){}

		counting_listener(Listener const &s)
		:Listener(s)
		,numberOfTests(0),successfulTests(0),failedTests(0),errors(0),numberOfSuites(0){}

		void begin(suite const &s, char const *info){
			++numberOfSuites;
			Listener::begin(s,info);
		}
		void start(test const &t){
			++numberOfTests;
			Listener::start(t);
		}
		void success(test const &t,char const *msg){
			++successfulTests;
			Listener::success(t,msg);
		}
		void failure(test const &t,test_failure const &e){
			++failedTests;
			Listener::failure(t,e);
		}
		void error(test const &t,char const *what){
			++errors;
			Listener::error(t,what);
		}
		int numberOfTests;
		int successfulTests;
		int failedTests;
		int errors;
		int numberOfSuites;
	};
}
#endif /*CUTE_COUNTING_LISTENER_H_*/
