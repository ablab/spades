#include "logger/log_writers.hpp"

#undef INFO
#define INFO(message)                       \
{                                                                         \
    cout << __FILE__ << " " <<  __LINE__ << "  :::  " << message << endl; \
}                                                                         \


#define LOG(message)                                                      \
{                                                                         \
    cout << __FILE__ << " " <<  __LINE__ << "  :::  " << message << endl; \
}                                                                         \

