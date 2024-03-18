//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __CONFIG_SINGL_HPP__
#define __CONFIG_SINGL_HPP__

#include "utils/verify.hpp"

#include <string>

namespace config_common {

// config singleton-wrap
template<class Config>
struct config {

    template<class Source>
    static void create_instance(Source const &source) {
        load(inner_cfg(), source);
        is_initialized() = true;
    }

    static Config const &get() {
        VERIFY_MSG(is_initialized(), "Config not initialized");
        return inner_cfg();
    }

    static Config &get_writable() {
        VERIFY_MSG(is_initialized(), "Config not initialized");
        return inner_cfg();
    }

private:
    static Config &inner_cfg() {
        static Config config;
        return config;
    }

    static bool &is_initialized() {
        static bool is_initialized = false;
        return is_initialized;
    }
};

}


#endif // __CONFIG_SINGLE_HPP__
