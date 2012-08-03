//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#if __DARWIN || __DARWIN_UNIX03
#include <mach/task.h>
#include <mach/mach.h>
#else
#include <sys/resource.h>
#endif

namespace logging
{

#if __DARWIN || __DARWIN_UNIX03
inline unsigned get_max_rss() {
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if (KERN_SUCCESS !=
      task_info(mach_task_self(),
                TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count))
    return -1U;

  return t_info.resident_size / 1024;
}
#else
inline unsigned get_max_rss() {
  rusage ru;
  getrusage(RUSAGE_SELF, &ru);

  return ru.ru_maxrss;
}
#endif

inline properties::properties(level default_level)
    : def_level(default_level)
{
}

inline properties::properties(string filename, level default_level)
    : def_level(default_level)
{
    if (filename.empty())
        return;

    ifstream in(filename);

    map<string, level> remap =
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

        string str(buf);
        trim(str);

        if (str.empty() || boost::starts_with(str, "#"))
            continue;

        vector<string> entry;
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

    auto it = props_.levels.find(source);
    if (it != props_.levels.end())
        source_level = it->second;

    return desired_level >= source_level;
}

inline void logger::log(level desired_level, const char* file, size_t line_num, const char* source, const char* msg)
{
    double time = timer_.time();
    unsigned max_rss = get_max_rss();

    for(auto it = writers_.begin(); it != writers_.end(); ++it)
        (*it)->write_msg(time, max_rss, desired_level, file, line_num, source, msg);
}

//
inline void logger::add_writer(writer_ptr ptr)
{
    writers_.push_back(ptr);
}

////////////////////////////////////////////////////
inline optional<logger>& __logger()
{
    static optional<logger> l;
    return l;
}

inline void create_logger(string filename, level default_level)
{
    properties props(filename, default_level);
    __logger() = in_place(props);
}

} // logging
