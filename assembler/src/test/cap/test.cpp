// just to check that headers from include and debruijn folders are correctly included
#include "standard.hpp"
#include "logger/log_writers.hpp"
#include "segfault_handler.hpp"
#include "stacktrace.hpp"
#include "config_struct.hpp"
#include "io/easy_reader.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/cutting_reader_wrapper.hpp"
#include "io/multifile_reader.hpp"
#include "io/careful_filtering_reader_wrapper.hpp"
#include "launch.hpp"
#include "simple_tools.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "omni/distance_estimation.hpp"
#include "memory_limit.hpp"
#include "boost/archive/tmpdir.hpp"
#include "read_converter.hpp"

int main() {
    std::cout << "Hello world! I'm in CAP TEST project" << std::endl;
    return 0;
}
