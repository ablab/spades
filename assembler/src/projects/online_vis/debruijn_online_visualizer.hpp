//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "online_visualizer.hpp"
#include "debruijn_environment.hpp"
#include "debruijn_commands.hpp"

namespace online_visualization {

class DebruijnOnlineVisualizer : public OnlineVisualizer<DebruijnEnvironment> {
 protected:
  void AddSpecificCommands() {
    AddCommand(make_shared<LoadGenomeCommand>());
    AddCommand(make_shared<SetMaxVertCommand>());
    AddCommand(make_shared<SetFolderCommand>());
    AddCommand(make_shared<SetFileNameCommand>());

    AddCommand(make_shared<FillPositionCommand>());
    AddCommand(make_shared<ClearPositionCommand>());

    AddCommand(make_shared<DrawVertexCommand>());
    AddCommand(make_shared<DrawEdgeCommand>());
    AddCommand(make_shared<DrawPositionCommand>());
    AddCommand(make_shared<DrawPartOfGenomeCommand>());
    AddCommand(make_shared<DrawContigCommand>());
    AddCommand(make_shared<DrawContigsCommand>());
    AddCommand(make_shared<DrawPolymorphicRegions>());
    AddCommand(make_shared<DrawPoorlyAssembledCommand>());
    AddCommand(make_shared<DrawUnresolvedWRTAssemblyCommand>());
    AddCommand(make_shared<DrawUnresolvedWRTReferenceCommand>());
    AddCommand(make_shared<DrawConnectedCommand>());
    AddCommand(make_shared<ShowPositionCommand>());
    AddCommand(make_shared<DrawMisassemblies>());
    AddCommand(make_shared<DrawCoverageDropsCommand>());
    AddCommand(make_shared<PrintPathsCommand>());
    AddCommand(make_shared<PrintContigsStatsCommand>());
    AddCommand(make_shared<JunctionSequenceCommand>());
    AddCommand(make_shared<PrintEdgeCommand>());
    AddCommand(make_shared<ClipTipsCommand>());
  }

 public:
  DebruijnOnlineVisualizer() : OnlineVisualizer<DebruijnEnvironment>() {
  }
};

}
