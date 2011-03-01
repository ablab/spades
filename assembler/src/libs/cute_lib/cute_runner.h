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

#ifndef CUTE_RUNNER_H_
#define CUTE_RUNNER_H_
#include "cute_test.h"
#include "cute_suite.h"
#include "cute_listener.h"
namespace cute {
	template <typename Listener=null_listener>
	struct runner : Listener{
		runner():Listener(){}
		runner(Listener &s):Listener(s){}
		bool operator()(test const &t){
			return runit(t);
		}
		bool operator()(suite const &s,char const *info=""){
			Listener::begin(s,info);
			bool result=true;
			for(suite::const_iterator it=s.begin();
			    it != s.end();
			    ++it){
			    	result = this->runit(*it) && result;
			    }
			Listener::end(s,info);
			return result;
		}
	private:
		bool runit(test const &t){
			try {
				Listener::start(t);
				t();
				Listener::success(t,"OK");
				return true;
			} catch (cute::test_failure const &e){
				Listener::failure(t,e);
			} catch (std::exception const &exc){
				Listener::error(t,demangle(exc.what()).c_str());
			} catch (std::string &s){
				Listener::error(t,s.c_str());
			} catch (char const *&cs) {
				Listener::error(t,cs);
			} catch(...) {
				Listener::error(t,"unknown exception thrown");
			}
			return false;
		}
	};
	template <typename Listener>
	runner<Listener> makeRunner(Listener &s){
		return runner<Listener>(s);
	}
}
#endif /*CUTE_RUNNER_H_*/
