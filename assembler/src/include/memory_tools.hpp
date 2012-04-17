//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <sys/time.h>
#include <sys/resource.h>

bool print_mem_usage(std::string const& msg)
{
	static size_t pid = getpid();
	string str = (format("pmap -d %d | grep writeable/private") % pid).str();
	cout << "==== MEM USAGE: " << msg << endl;
	return system(str.c_str()) == 0;
}

void limit_memory(size_t limit){

	rlimit rl = {limit, limit};

	if (sizeof(rl.rlim_max) < 8){

		INFO("Can't limit virtual memory because of 32-bit system");
		return;
	}

	int res = setrlimit(RLIMIT_AS, &rl);
    VERIFY(res == 0);
}
