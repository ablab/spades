//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

namespace logging
{

inline properties::properties()
    : def_level(L_INFO)
{
}

inline properties::properties(string filename)
    : def_level(L_INFO)
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

    for(auto it = writers_.begin(); it != writers_.end(); ++it)
        (*it)->write_msg(time, desired_level, file, line_num, source, msg);
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

inline void create_logger(string filename = "")
{
    properties props(filename);
    __logger() = in_place(props);
}

} // logging
