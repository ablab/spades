//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef HAMMER_STATS_HPP
#define HAMMER_STATS_HPP

namespace hammer {
void print_time();
void print_mem_usage();
void print_stats();
void print_full_stats();
};

#define TIMEDLN(a)                                            \
  do {                                                        \
    hammer::print_full_stats(); std::cout << a << std::endl;  \
  } while (0)

#endif // HAMMER_STATS_HPP
