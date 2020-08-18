//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "decl_logger.hpp"
#include "utils/perf/perfcounter.hpp"
#include "version.hpp"

#include <vector>
#include <unordered_map>
#include <string>
#include <sstream>
#include <memory>
#include <fstream>
#include <iostream>

#include "config.hpp"

namespace logging
{

/////////////////////////////////////////////////////
enum level
{
    L_TRACE,
    L_DEBUG,
    L_INFO,
    L_WARN,
    L_ERROR
};

inline std::string level_name(level l)
{
  static std::string names [] =
    {
        "TRACE",
        "DEBUG",
        "INFO" ,
        "WARN" ,
        "ERROR"
    };

    return names[l];
}


/////////////////////////////////////////////////////
struct writer
{
  virtual void write_msg(double time_in_sec, size_t cmem, size_t max_rss, level l, const char* file, size_t line_num, const char* source, const char* msg) = 0;

  virtual ~writer(){}
};

typedef std::shared_ptr<writer> writer_ptr;

/////////////////////////////////////////////////////
struct properties
{
    /* Reading logger properties from file
     *
     * File should contains lines like below.
     * Use leading # for comment.
     * File could contain line with default behavior description. If no 'default' entry found, default is set to INFO
     * Valid levels: TRACE, DEBUG, INFO, WARN, ERROR
     *
     *    default=INFO
     *    AbraCaDabra=TRACE
     *    #BubaZuba=WARN
     *    HariKrishna=INFO
     *
     */

    properties(std::string filename = "", level default_level = L_INFO);
    properties(level default_level = L_INFO);

    std::unordered_map<std::string, level> levels;
    level  def_level;
    bool   all_default;
};

////////////////////////////////////////////////////
struct logger
{
    logger(properties const& props);

    //
    bool need_log(level desired_level, const char* source) const;
    void log(level desired_level, const char* file, size_t line_num, const char* source, const char* msg);

    //
    void add_writer(writer_ptr ptr) {
        writers_.push_back(ptr);
    }

    template<class Writer, typename... Args>
    void add_writer(Args&&... args) {
        writers_.push_back(std::make_shared<Writer>(std::forward<Args>(args)...));
    }
    

private:
    properties                 props_  ;
    std::vector<writer_ptr>    writers_;
    utils::perf_counter        timer_  ;
};

std::shared_ptr<logger>& __logger();
logger* create_logger(std::string filename = "", level default_level = L_INFO);

void attach_logger(logger *lg);
void detach_logger();

} // logging

#define LOG_MSG(l, msg)                                                 \
  do {                                                                  \
    std::shared_ptr<logging::logger> &__lg__ = logging::__logger();     \
    if (__lg__.get() == NULL) {                                         \
      std::cerr << "WARNING: Try to use logger before create one. Level=" << level_name(l);\
      std::cerr << ". Message="  <<  msg << "\n";                       \
      fflush(stderr);                                                   \
      break;                                                            \
    }                                                                   \
                                                                        \
    if (__lg__->need_log((l), __scope_source_name())) {                 \
      std::stringstream __logger__str__;                                \
      __logger__str__ << msg; /* don't use brackets here! */            \
      __lg__->log((l), __FILE__, __LINE__, __scope_source_name(), __logger__str__.str().c_str()); \
    }                                                                   \
  } while(0);

#define LOG_EXPR(l, expr)                                               \
    do {                                                                \
        std::shared_ptr<logging::logger> &__lg__ = logging::__logger(); \
        if (__lg__.get() == NULL) {                                     \
            std::cerr << "WARNING: Try to use logger before create one. Level=" << level_name(l); \
            fflush(stderr);                                             \
            break;                                                      \
        }                                                               \
        if (__lg__->need_log((l), __scope_source_name())) {             \
            expr;                                                       \
        }                                                               \
    } while (0);

#ifdef SPADES_DEBUG_LOGGING
# define DEBUG(message)                     LOG_MSG(logging::L_DEBUG, message)
# define DEBUG_EXPR(expr)                   LOG_EXPR(logging::L_DEBUG, expr)
# define TRACE(message)                     LOG_MSG(logging::L_TRACE, message)
# define TRACE_EXPR(expr)                   LOG_EXPR(logging::L_TRACE, expr)
#else
# define DEBUG(message)                     /* No trace */
# define DEBUG_EXPR(expr)                   /* No trace */
# define TRACE(message)                     /* No trace */
# define TRACE_EXPR(expr)                   /* No trace */
#endif
#define INFO(message)                       LOG_MSG(logging::L_INFO , message)
#define START_BANNER(description)                                       \
    do {                                                                \
        INFO("Starting " description ", built from "                    \
             << version::refspec()                                      \
             << ", git revision "                                       \
             << version::gitrev());                                     \
    } while (0)
#define VERBOSE_T(n, T, message)            {size_t n_copy = (n); if (n_copy % (T) == 0 && n_copy > 0) INFO(n_copy << message)}
#define VERBOSE(n, message)                 VERBOSE_T((n), 10000, message)
#define VERBOSE_POWER_T(n, T, message)      {size_t n_copy = (n); if ((n_copy & (n_copy - 1)) == 0 && (n_copy > T)) INFO(n_copy << message)}
#define VERBOSE_POWER(n, message)           VERBOSE_POWER_T((n), 10000, message)
#define VERBOSE_POWER_T2(n, T, message)     {size_t n_copy = (n); if ((n_copy & (n_copy - 1)) == 0 && (n_copy > T)) INFO(message)}
#define VERBOSE_POWER2(n, message)          VERBOSE_POWER_T2((n), 10000, message)
#define WARN(message)                       LOG_MSG(logging::L_WARN, message)
#define ERROR(message)                      LOG_MSG(logging::L_ERROR, message)
#define FATAL_ERROR(message)                                            \
    do {                                                                \
        ERROR(message);                                                 \
        if (errno != 0) {                                               \
            exit(errno);                                                \
        } else {                                                        \
            exit(-1);                                                   \
        }                                                               \
    } while(0)

#define CHECK_FATAL_ERROR(expr, msg)                                    \
    if (!(expr)) {                                                      \
        FATAL_ERROR(msg);                                               \
    }
