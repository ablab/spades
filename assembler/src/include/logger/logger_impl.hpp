//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include <string>
#include <map>
#include <fstream>
#include <vector>

#include "config.hpp"

#ifdef SPADES_USE_JEMALLOC
# include <jemalloc/jemalloc.h>
#endif

namespace logging
{

inline properties::properties(level default_level)
    : def_level(default_level), all_default(true)
{
}

inline properties::properties(std::string filename, level default_level)
    : def_level(default_level), all_default(true)
{
    if (filename.empty())
        return;

    std::ifstream in(filename.c_str());

    std::map<std::string, level> remap =
    {
        {"TRACE", L_TRACE},
        {"DEBUG", L_DEBUG},
        {"INFO" , L_INFO },
        {"WARN" , L_WARN },
        {"ERROR", L_ERROR}
    };

    while (!in.eof())
    {
        using namespace boost;

        char buf [0x400] = {};
        in.getline(buf, sizeof buf);

        std::string str(buf);
        trim(str);

        if (str.empty() || boost::starts_with(str, "#"))
            continue;

        std::vector<std::string> entry;
        split(entry, str, is_any_of("="));

        if(entry.size() != 2)
            throw std::runtime_error("invalid log file property entry: " + str);

        trim    (entry[0]);
        trim    (entry[1]);
        to_upper(entry[1]);

        auto it = remap.find(entry[1]);
        if(it == remap.end())
            throw std::runtime_error("invalid log file level description: " + entry[1]);

        levels[entry[0]] = it->second;
    }

    auto def = levels.find("default");
    if (def != levels.end())
        def_level = def->second;

    for (auto I = levels.begin(), E = levels.end(); I != E; ++I) {
      if (I->second != def_level) {
        all_default = false;
        break;
      }
    }
}


////////////////////////////////////////////////////
inline logger::logger(properties const& props)
    : props_(props)
{
}

//
inline bool logger::need_log(level desired_level, const char* source) const
{
    level source_level = props_.def_level;

    if (!props_.all_default) {
      auto it = props_.levels.find(source);
      if (it != props_.levels.end())
        source_level = it->second;
    }

    return desired_level >= source_level;
}

#ifdef SPADES_USE_JEMALLOC

inline void logger::log(level desired_level, const char* file, size_t line_num, const char* source, const char* msg) {
  double time = timer_.time();
  const size_t *cmem = 0, *cmem_max = 0;
  size_t clen = sizeof(cmem);

  je_mallctl("stats.cactive", &cmem, &clen, NULL, 0);
  je_mallctl("stats.cactive_max", &cmem_max, &clen, NULL, 0);

  for (auto it = writers_.begin(); it != writers_.end(); ++it)
    (*it)->write_msg(time, (*cmem) / 1024, (*cmem_max) / 1024, desired_level, file, line_num, source, msg);
}
#else
inline void logger::log(level desired_level, const char* file, size_t line_num, const char* source, const char* msg) {
  double time = timer_.time();
  size_t max_rss = get_max_rss();

  for (auto it = writers_.begin(); it != writers_.end(); ++it)
    (*it)->write_msg(time, max_rss, desired_level, file, line_num, source, msg);
}
#endif

//
inline void logger::add_writer(writer_ptr ptr)
{
    writers_.push_back(ptr);
}

////////////////////////////////////////////////////
inline std::shared_ptr<logger> &__logger() {
  static std::shared_ptr<logger> l;
  return l;
}

inline logger *create_logger(std::string filename, level default_level) {
  return new logger(properties(filename, default_level));
}

inline void attach_logger(logger *lg) {
  __logger().reset(lg);
}

inline void detach_logger() {
  __logger().reset();
}


} // logging
