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


    Example of extending a type.
    Usage:
    	short option: -d mm/dd/yyyy
    	long option:  --date mm/dd/yyyy
    
*/

#include <iostream>
#include "getopt_pp.h"

using namespace GetOpt;

struct Date
{
	unsigned int year;
	unsigned int month;
	unsigned int day;
	
	bool valid() const { return (month >= 1 && month <= 12 && day >= 1 && day <= 31); }
	Date(){}
	Date(unsigned int y, unsigned int m, unsigned int d) : year(y), month(m), day(d) {}
		

};

namespace GetOpt
{
	template <> _Option::Result convert<Date>(const std::string& s, Date& d, std::ios::fmtflags)
	{
		_Option::Result ret = _Option::BadType;
		Date tmp;
		char slash;
		std::stringstream ss(s);
		if (ss >> tmp.month)
		{
			ss >> slash;
			if (ss >> tmp.day)
			{
				ss >> slash;
				if (ss >> tmp.year)
				{
					if (tmp.valid())
					{
						ret = _Option::OK;
						d = tmp;
					}
				}
			}
		}
		
		return ret;
	}
}

int main(int argc, char* argv[])
{
	Date date;
	const Date myBirthday(1977, 7, 31);
	
	GetOpt_pp ops(argc, argv);
	
	if (ops >> Option('d', "date", date, myBirthday))
		std::cout << "Date: " << date.month << "-" << date.day << "-" << date.year << std::endl;
	else
		std::cerr << "Invalid date. Enter mm/dd/yyyy" << std::endl;
	
	return 0;
}

