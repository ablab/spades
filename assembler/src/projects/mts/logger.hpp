#include "utils/logger/log_writers.hpp"

void create_console_logger() {
    using namespace logging;

    std::filesystem::path log_props_file = "log.properties";

    logger *lg = create_logger(exists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}
