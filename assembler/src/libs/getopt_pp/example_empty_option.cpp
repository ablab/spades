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
    

    Example testing empty options.
    Usage:
    	pass a list of arguments to the program, without specifying any option.
    	For example: ./program one two three
    
*/

#include <iostream>
#include "getopt_pp.h"

using namespace GetOpt;

int main(int argc, char* argv[])
{
	std::vector<std::string> args;
	
	GetOpt_pp ops(argc, argv);
	
	ops >> Option(GetOpt_pp::EMPTY_OPTION,args );
	
	std::cout << "RAN: " << ops.app_name() << " ";
	for (std::vector<std::string>::const_iterator it = args.begin(); it != args.end(); ++it)
		std::cout << *it << " ";
		
	std::cout<<std::endl;
	
	return 0;
}

