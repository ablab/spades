//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../online_vis/online_visualizer.hpp"
#include "cap_environment.hpp"
#include "cap_commands.hpp"

namespace online_visualization {

class CapOnlineVisualizer : public OnlineVisualizer<CapEnvironment> {
 protected:
  void AddSpecificCommands() {
    AddCommand(make_shared<AddGenomeCommand>());
    AddCommand(make_shared<BuildGraphCommand>());
    AddCommand(make_shared<RefineCommand>());
    AddCommand(make_shared<SaveGenomesCommand>());
    AddCommand(make_shared<SaveGraphCommand>());
    AddCommand(make_shared<SaveEnvCommand>());
    AddCommand(make_shared<LoadEnvCommand>());
    AddCommand(make_shared<DrawPicsCommand>());
    AddCommand(make_shared<FindIndelsCommand>());
    AddCommand(make_shared<FindInversionsCommand>());
    AddCommand(make_shared<SaveBlocksCommand>());
    AddCommand(make_shared<MosaicAnalysisCommand>());
    AddCommand(make_shared<MaskRepeatsCommand>());
    AddCommand(make_shared<LoadGraphCommand>());
    AddCommand(make_shared<BlocksToGRIMMFormat>());
  }

 public:
  CapOnlineVisualizer() : OnlineVisualizer<CapEnvironment>() {
  }
};

}
