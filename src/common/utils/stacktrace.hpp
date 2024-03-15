//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * stacktrace.hpp
 *
 *  Created on: Feb 2, 2012
 *      Author: valery
 */

#pragma once
#include <execinfo.h>
#include <iostream>

namespace utils {

inline void print_stacktrace() {
    std::cout << "=== Stack Trace ===" << std::endl;

    const size_t max_stack_size = 1000;

    void *stack_pointers[max_stack_size];
    int count = backtrace(stack_pointers, max_stack_size);

    char **func_names = backtrace_symbols(stack_pointers, count);

    // Print the stack trace
    for (int i = 0; i < count; ++i)
        std::cerr << func_names[i] << std::endl;

    // Free the string pointers
    free(func_names);
}

}
