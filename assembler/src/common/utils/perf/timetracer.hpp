//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#if 1
#include <llvm/Support/TimeProfiler.h>

namespace {
struct time_trace_scope {
    time_trace_scope()
            : trace(__PRETTY_FUNCTION__) {}

    time_trace_scope(llvm::StringRef comment)
            : trace(comment) {}

    time_trace_scope(llvm::StringRef comment, llvm::StringRef detail)
            : trace(comment, detail) {}

    llvm::TimeTraceScope trace;
};
};

#define _GET_OVERRIDE(_1, _2, _3, NAME, ...) NAME

#define TIME_TRACE_SCOPE_IMPL1(suf)  time_trace_scope trace ## suf
#define TIME_TRACE_SCOPE_IMPL2(suf, comment)  time_trace_scope trace ## suf(comment)
#define TIME_TRACE_SCOPE_IMPL3(suf, comment, detail)  time_trace_scope trace ## suf(comment, detail)
#define TIME_TRACE_SCOPE_IMPL(...)              \
    _GET_OVERRIDE(__VA_ARGS__,                  \
                  TIME_TRACE_SCOPE_IMPL3, TIME_TRACE_SCOPE_IMPL2, TIME_TRACE_SCOPE_IMPL1)(__VA_ARGS__)
#define TIME_TRACE_SCOPE(...)  TIME_TRACE_SCOPE_IMPL(__LINE__, ##__VA_ARGS__)
#define TIME_TRACE_BEGIN(comment) do {                              \
        llvm::timeTraceProfilerBegin(comment, llvm::StringRef("")); \
    } while(0);
#define TIME_TRACE_END do {                                         \
        llvm::timeTraceProfilerEnd();                               \
    } while(0);
    
#else
#define TIME_TRACE_SCOPE
#define TIME_TRACE_BEGIN
#define TIME_TRACE_END
#endif

