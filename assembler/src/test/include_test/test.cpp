//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
#include "seq_test.hpp"
#include "rtseq_test.hpp"
#include "sequence_test.hpp"
#include "quality_test.hpp"
#include "nucl_test.hpp"
#include "cyclic_hash_test.hpp"
#include "binary_test.hpp"*/

#include "utils/segfault_handler.hpp"
#include "utils/logger/logger.hpp"
#include "utils/logger/log_writers.hpp"

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
