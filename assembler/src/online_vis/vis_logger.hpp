#include "logger/log_writers.hpp"

#undef INFO
#define INFO(message)                       \
{                                                                         \
    cout << __FILE__ << " " <<  __LINE__ << "  :::  " << message << endl; \
}                                                                         \


#define LOG(message)                                                      \
{                                                                         \
    cout << message << endl;                                              \
}                                                                         \

//#define trace(message)                      LOG_MSG(logging::L_TRACE, message)
#define debug(print, message)               \
{                                           \
    if (print) {                            \
        cout << message << endl;            \
    }                                       \
}                                           
