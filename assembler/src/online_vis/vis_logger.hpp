#include "logger/log_writers.hpp"

#undef INFO
#define INFO(message)                       LOG_MSG(logging::L_DEBUG , message)

#define LOG(message)                                                      \
{                                                                         \
    cout << __FILE__ << " " <<  __LINE__ << "  :::  " << message << endl; \
}                                                                         \

