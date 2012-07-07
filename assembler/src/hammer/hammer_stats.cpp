//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <boost/format.hpp>
#include "hammer_stats.hpp"

namespace hammer {

#if __DARWIN || __DARWIN_UNIX03
#include <mach/task.h>
#include <mach/mach.h>

unsigned get_max_rss() {
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if (KERN_SUCCESS !=
      task_info(mach_task_self(),
                TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count))
    return -1U;

  return t_info.resident_size / 1024;
}
#else
#include <sys/resource.h>

unsigned get_max_rss() {
  rusage ru;
  getrusage(RUSAGE_SELF, &ru);

  return ru.ru_maxrss;
}
#endif

void print_time() {
  time_t rawtime;
  tm *ptm;
  time(&rawtime);
  ptm = gmtime(&rawtime);
  std::cout << std::setfill('0') << "[ " << std::setw(2) << ptm->tm_hour << ":" << std::setw(2) << ptm->tm_min << ":" << std::setw(2) << ptm->tm_sec << " ] ";
}

void print_mem_usage() {
  static size_t pid = getpid();
  std::string str = (boost::format("pmap -d %d | grep writeable/private") % pid).str();
  //std::cout << "==== MEM USAGE ==== " << std::std::endl;
  if (system(str.c_str()))
    std::cout << "  System error!" << std::endl;
}

void print_stats() {
  std::cout << "[";
  time_t rawtime;
  tm * ptm2;
  time(&rawtime);
  ptm2 = gmtime(&rawtime);
  std::cout << " " << std::setfill('0') << std::setw(2) << ptm2->tm_hour << ":" << std::setw(2) << ptm2->tm_min << ":" << std::setw(2) << ptm2->tm_sec;

  rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  tm * ptm = gmtime(&ru.ru_utime.tv_sec);
  unsigned max_rss = get_max_rss();
  std::cout << " " << std::setw(2) << std::setfill('0')  << ptm->tm_hour << ":" << std::setw(2) << std::setfill('0') << ptm->tm_min << ":" << std::setw(2) << std::setfill('0') << ptm->tm_sec;
  if (max_rss < 1024 * 1024) {
    std::cout << std::setw(5) << std::setfill(' ') << (max_rss / 1024) << "M ";
  } else {
    std::cout << std::setw(6) << std::setprecision(1) << std::fixed << std::setfill(' ') << (max_rss / (1024.0 * 1024.0) ) << "G ";
  }
  std::cout << "] ";
}

void print_full_stats() {
  // print_mem_usage();
  std::cout << "[";
  time_t rawtime;
  tm *ptm2;
  time (&rawtime);
  ptm2 = gmtime(&rawtime);
  std::cout << " " << std::setfill('0') << std::setw(2) << ptm2->tm_hour << ":" << std::setw(2) << ptm2->tm_min << ":" << std::setw(2) << ptm2->tm_sec;

  rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  tm * ptm = gmtime(&ru.ru_utime.tv_sec);
  unsigned max_rss = get_max_rss();
  std::cout << " " << std::setw(2) << std::setfill('0')  << ptm->tm_hour << ":" << std::setw(2) << std::setfill('0') << ptm->tm_min << ":" << std::setw(2) << std::setfill('0') << ptm->tm_sec;
  if (max_rss < 1024 * 1024) {
    std::cout << std::setw(5) << std::setfill(' ') << (max_rss / 1024) << "M ";
  } else {
    std::cout << std::setw(6) << std::setprecision(1) << std::fixed << std::setfill(' ') << (max_rss / (1024.0 * 1024.0) ) << "G ";
  }
  std::cout << "] ";
}

};
