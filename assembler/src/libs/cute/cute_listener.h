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

#ifndef CUTE_LISTENER_H_
#define CUTE_LISTENER_H_
#include "cute_base.h"
#include "cute_test.h"
#include "cute_suite.h"
namespace cute {
	struct null_listener{ // defines Contract of runner parameter
		void begin(suite const &s, char const *info){}
		void end(suite const &s, char const *info){}
		void start(test const &t){}
		void success(test const &t,char const *msg){}
		void failure(test const &t,test_failure const &e){}
		void error(test const &t,char const *what){}
	};
}
#endif /* CUTE_LISTENER_H_ */

