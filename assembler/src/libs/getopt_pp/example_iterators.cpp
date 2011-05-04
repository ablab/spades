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


    Example of iterating over the options.
    Usage:
    	nothing special. Just invoke this with any options, it will dump them.
*/

#include <iostream>
#include "getopt_pp.h"

using namespace GetOpt;

int main(int argc, char* argv[])
{
	GetOpt_pp ops(argc, argv);
	char short_opt[2] = {0, 0};

	std::cout << "Short options:" << std::endl;
	
	for(GetOpt_pp::short_iterator it = ops.begin(); it != ops.end(); ++it)
	{
		short_opt[0] = it.option();
		std::cout << "\t" << short_opt << " has " << it.args().size() << " arguments." << std::endl;
	}	
	
	std::cout << std::endl << "Long options:" << std::endl;
	
	for(GetOpt_pp::long_iterator it = ops.begin(); it != ops.end(); ++it)
		std::cout << "\t" << it.option() << " has " << it.args().size() << " arguments." << std::endl;
	
	return 0;
}

