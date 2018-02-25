#include <gtest/gtest.h>
#include <stdio.h>
#include "utils/segfault_handler.hpp"

GTEST_API_ int main(int argc, char **argv) {
  utils::segfault_handler sh;
  printf("Running main() from gtest_main.cpp\n");
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
