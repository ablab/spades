//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************


#pragma once

#include "utils/stacktrace.hpp"
#include "boost/noncopyable.hpp"

#include <functional>
#include <signal.h>

namespace utils {

struct segfault_handler : boost::noncopyable {
    typedef std::function<void()> callback_t;

    typedef void (*seg_handler_t)(int);

    segfault_handler(callback_t const &cb = 0) {
        if (callback() != 0)
            throw std::runtime_error("failed to initialize segfault_handler, it has been already initialized");

        callback() = cb;
        old_func_ = signal(SIGSEGV, &segfault_handler::handler);
    }

    ~segfault_handler() {
        callback() = 0;
        signal(SIGSEGV, old_func_);
    }

private:
    static callback_t &callback() {
        static callback_t cb = 0;
        return cb;
    }

    static void handler(int signum) {
        if (signum == SIGSEGV) {
            std::cerr << "The program was terminated by segmentation fault" << std::endl;
            print_stacktrace();

            if (callback())
                callback()();
        }

        //TEST!!
        exit(1);

        signal(signum, SIG_DFL);
        kill(getpid(), signum);
    }

private:
    seg_handler_t old_func_;
};

}
