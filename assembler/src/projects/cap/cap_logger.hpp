//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/logger/log_writers.hpp"

/*
#undef INFO
#define INFO(message)                       \
{                                                                         \
    cout << __FILE__ << " " <<  __LINE__ << "  :::  " << message << endl; \
}                                                                         \
*/


#define LOG(message)                                                      \
{                                                                         \
    cout << message << endl;                                              \
}                                                                         \

//#define trace(message)                      LOG_MSG(logging::L_TRACE, message)
#define debug(print, message)               \
{                                           \
    if (print) {                            \
        cout << message << endl;            \
    }                                       \
}                                           
