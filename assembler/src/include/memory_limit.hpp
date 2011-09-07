#pragma once

#include <sys/time.h>
#include <sys/resource.h>

void limit_memory(size_t limit){

    rlimit rl = {limit, limit};
    int res = setrlimit(RLIMIT_AS, &rl);

    assert(res == 0);
}
