#include <gtest/gtest.h>
#include <stdio.h>
#include "utils/segfault_handler.hpp"
#include "utils/logger/logger.hpp"
#include "utils/logger/log_writers.hpp"
#include "utils.hpp"

GTEST_API_ int main(int argc, char **argv) {
  utils::segfault_handler sh;
  create_console_logger();
  printf("Running main() from gtest_main.cpp\n");
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
