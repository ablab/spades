//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/stacktrace.hpp"

#include <sstream>
#include <iostream>
#include <cassert>

#include "config.hpp"

#define VERIFY(expr)                                             \
    do {                                                         \
        if (!(expr))                                             \
            utils::print_stacktrace();                           \
        assert(expr);                                            \
    } while(0);

#ifdef SPADES_ENABLE_EXPENSIVE_CHECKS
# define VERIFY_DEV(expr) VERIFY(expr)
#else
# define VERIFY_DEV(expr) (void)(true ? (void)0 : ((void)(expr)));
#endif

#define VERIFY_MSG(expr, msg)                                           \
    if (!(expr)) {                                                      \
        std::stringstream ss;                                           \
        utils::print_stacktrace();                                      \
        ss << "Verification of expression '" << #expr << "' failed in function '" <<  __PRETTY_FUNCTION__ << \
                "'. In file '" << __FILE__ << "' on line " << __LINE__ << ". Message '" << msg << "'." ; \
        std::cout << ss.str() << std::endl;                             \
        std::cerr << ss.str() << std::endl;                             \
        fflush(stdout);                                                 \
        fflush(stderr);                                                 \
        assert(expr);                                                   \
    }
