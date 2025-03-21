//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "config.hpp"

#include "utils/logger/logger.hpp"
#include "utils/memory_limit.hpp"
#include "utils/parallel/openmp_wrapper.h"
#include "utils/perf/memory.hpp"
#include "utils/stl_utils.hpp"

#include <cppformat/format.h>

#include <fstream>
#include <map>
#include <string>
#include <vector>

#ifdef SPADES_USE_JEMALLOC
# include <jemalloc/jemalloc.h>
#endif

namespace logging {

properties::properties(level default_level)
        : def_level(default_level), all_default(true) {}

properties::properties(std::filesystem::path filename, level default_level)
    : def_level(default_level), all_default(true) {
    if (filename.empty())
        return;

    std::ifstream in(filename);

    std::map<std::string, level> remap = {
        {"TRACE", L_TRACE},
        {"DEBUG", L_DEBUG},
        {"INFO" , L_INFO },
        {"WARN" , L_WARN },
        {"ERROR", L_ERROR}
    };

    while (!in.eof()) {
        char buf [0x400] = {};
        in.getline(buf, sizeof buf);

        std::string str(buf);
        utils::trim(str);

        if (str.empty() || utils::starts_with(str, "#"))
            continue;

        auto res = utils::split(str, "=");
        if (res.size() != 2)
            throw std::runtime_error("invalid log file property entry: " + str);

        std::vector<std::string> entry(res.begin(), res.end());
        utils::trim(entry[0]);
        utils::trim(entry[1]);
        entry[1] = utils::str_toupper(entry[1]);

        auto it = remap.find(entry[1]);
        if (it == remap.end())
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


logger::logger(properties const& props)
    : props_(props) { }

bool logger::need_log(level desired_level, const char* source) const {
    level source_level = props_.def_level;

    if (!props_.all_default) {
      auto it = props_.levels.find(source);
      if (it != props_.levels.end())
        source_level = it->second;
    }

    return desired_level >= source_level;
}

void logger::log(level desired_level, const std::filesystem::path& file, size_t line_num, const char* source, const char* msg) {
  double time = timer_.time();
  size_t mem = -1ull;
  size_t max_rss = -1ull;

#if defined(SPADES_USE_JEMALLOC)
  // Cannot use FATAL_ERROR here, we're inside logger

  // Update statisitcs cached by mallctl
  {
      uint64_t epoch = 1;
      size_t sz = sizeof(epoch);
      if (je_mallctl("epoch", &epoch, &sz, &epoch, sz) != 0) {
          fprintf(stderr, "mallctl() call failed, errno = %d", errno);
          exit(errno);
      }
  }

  {
      size_t cmem = 0;
      size_t clen = sizeof(cmem);

      int res = je_mallctl("stats.resident", &cmem, &clen, NULL, 0);
      if (res != 0) {
          fprintf(stderr, "mallctl() call failed, errno = %d", errno);
          exit(errno);
      }
      mem = (cmem + 1023)/ 1024;
  }
  max_rss = utils::get_max_rss();
  if (mem > max_rss)
      max_rss = mem;
#elif defined(SPADES_USE_MIMALLOC)
  mem = (utils::get_used_memory() + 1023) / 1024;
  // FIXME: we may need to refine here
  max_rss = utils::get_max_rss();
  if (mem > max_rss)
      max_rss = mem;
#else
  max_rss = utils::get_max_rss();
#endif

  for (auto it = writers_.begin(); it != writers_.end(); ++it)
    (*it)->write_msg(time, mem, max_rss, desired_level, file, line_num, source, msg);
}

////////////////////////////////////////////////////
std::shared_ptr<logger> &__logger() {
  static std::shared_ptr<logger> l;
  return l;
}

logger *create_logger(std::filesystem::path filename, level default_level) {
  return new logger(properties(filename, default_level));
}

void attach_logger(logger *lg) {
  __logger().reset(lg);
}

void detach_logger() {
  __logger().reset();
}


} // logging
