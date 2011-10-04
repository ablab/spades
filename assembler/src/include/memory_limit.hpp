#pragma once

#include <sys/time.h>
#include <sys/resource.h>

void limit_memory(size_t limit){

	rlimit rl = {limit, limit};

	if (sizeof(rl.rlim_max) < 8){

		INFO("Can't limit virtual memory because of 32-bit system");
		return;
	}

	int res = setrlimit(RLIMIT_AS, &rl);
    VERIFY(res == 0);
}
