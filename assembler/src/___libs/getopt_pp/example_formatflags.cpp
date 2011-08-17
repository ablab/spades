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


    Example using format flags: entering numbers in hex.
    Usage:
    	short option: -i number
    	long option:  --number number
    where number has the form 0xNNNN
    
*/

#include <iostream>
#include "getopt_pp.h"

using namespace GetOpt;

int main(int argc, char* argv[])
{
	int i;
	
	GetOpt_pp ops(argc, argv);
	
	ops >> std::hex >> Option('i', "number", i);
	
	std::cout << std::hex << "Number entered: (hex)" << i << " (dec)" << std::dec << i << std::endl;
	
	return 0;
}

