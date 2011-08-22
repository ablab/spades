/*
GetOpt_pp:	Yet another C++ version of getopt.
    Copyright (C) 2007, 2008  Daniel Gutson, FuDePAN
    
    This file is part of GetOpt_pp.

    GetOpt_pp is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    board-games is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    
	This file is a sample usage (and test case).
	TODO:
		* Support for arrays (of any type):
			opt >> vector<T>
		* More validations
		* Pick good exceptions types
		* Support the 'empty option' arguments
		* fill a structure at once
    
*/

#include <iostream>
#include "getopt_pp.h"

using namespace GetOpt;

int main(int argc, char* argv[])
{
	int test1 = 10;
	float test2 = 3.14f;
	std::string test3 = "hello";
	bool flag;
	std::vector<int> vec_int;
	
	try
	{	
		GetOpt_pp ops(argc, argv);
					
		ops.exceptions(std::ios::failbit | std::ios::eofbit);
		
		ops 
			>> Option('i', "test1", test1)
			>> Option('f', test2)
			>> Option('s', "string", test3)
			>> OptionPresent('x', "flag", flag)
			>> Option('v', "vec", vec_int)
		;
		
		if (!ops.options_remain())
		{
			std::cout << test1 << "\n" << test2 << "\n" << test3 << "\n" << flag << "\n";
			
			for(std::vector<int>::const_iterator it = vec_int.begin(); it != vec_int.end(); ++it)
				std::cout << *it << " ";
				
			std::cout << std::endl;

			return 0;
		}
		else
		{
			std::cerr << "too many options" << std::endl;
			return 1;
		}
	}
	catch(const GetOptEx& e)
	{
		std::cerr << "Invalid options\n";
	}
}

