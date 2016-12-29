//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/segfault_handler.hpp"
#include "utils/logger/logger.hpp"
#include "utils/logger/log_writers.hpp"
#include "utils/parallel/openmp_wrapper.h"

#include <gtest/gtest.h>
#include <teamcity_gtest/teamcity_gtest.h>
#include <cstdio>

void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

GTEST_API_ int main(int argc, char **argv) {
  utils::segfault_handler sh;
  create_console_logger();

  // Fix number of threads according to OMP capabilities.
  int max_threads = std::min(4, omp_get_max_threads());
  // Inform OpenMP runtime about this :)
  omp_set_num_threads(max_threads);

  printf("Running main() from gtest_main.cpp\n");
  testing::InitGoogleTest(&argc, argv);

  if (jetbrains::teamcity::underTeamcity()) {
      ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
      // Add unique flowId parameter if you want to run test processes in parallel
      // See http://confluence.jetbrains.net/display/TCD6/Build+Script+Interaction+with+TeamCity#BuildScriptInteractionwithTeamCity-MessageFlowId
      listeners.Append(new jetbrains::teamcity::TeamcityGoogleTestEventListener());
  }

  return RUN_ALL_TESTS();
}
