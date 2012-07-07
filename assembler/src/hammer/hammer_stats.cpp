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
  std::cout << " " << std::setw(2) << std::setfill('0')  << ptm->tm_hour << ":" << std::setw(2) << std::setfill('0') << ptm->tm_min << ":" << std::setw(2) << std::setfill('0') << ptm->tm_sec;
  if (ru.ru_maxrss < 1024 * 1024) {
    std::cout << std::setw(5) << std::setfill(' ') << (ru.ru_maxrss / 1024) << "M ";
  } else {
    std::cout << std::setw(6) << std::setprecision(1) << std::fixed << std::setfill(' ') << (ru.ru_maxrss / (1024.0 * 1024.0) ) << "G ";
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
  std::cout << " " << std::setw(2) << std::setfill('0')  << ptm->tm_hour << ":" << std::setw(2) << std::setfill('0') << ptm->tm_min << ":" << std::setw(2) << std::setfill('0') << ptm->tm_sec;
  if (ru.ru_maxrss < 1024 * 1024) {
    std::cout << std::setw(5) << std::setfill(' ') << (ru.ru_maxrss / 1024) << "M ";
  } else {
    std::cout << std::setw(6) << std::setprecision(1) << std::fixed << std::setfill(' ') << (ru.ru_maxrss / (1024.0 * 1024.0) ) << "G ";
  }
  std::cout << "] ";
}

};
