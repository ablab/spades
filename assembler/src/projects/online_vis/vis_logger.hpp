//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/logger/log_writers.hpp"

#undef INFO
#define INFO(message)                       \
{                                                                         \
    std::cout << __FILE__ << " " <<  __LINE__ << "  :::  " << message << std::endl; \
}                                                                         \


#define LOG(message)                                                      \
{                                                                         \
    std::cout << message << std::endl;                                              \
}                                                                         \

//#define trace(message)                      LOG_MSG(logging::L_TRACE, message)
#define debug(print, message)               \
{                                           \
    if (print) {                            \
        std::cout << message << std::endl;            \
    }                                       \
}                                           
