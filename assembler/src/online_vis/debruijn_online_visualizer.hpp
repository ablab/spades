#pragma once

#include "online_visualizer.hpp"
#include "debruijn_environment.hpp"
#include "debruijn_commands.hpp"

namespace online_visualization {

class DebruijnOnlineVisualizer : public OnlineVisualizer<DebruijnEnvironment> {
 protected:
  void AddSpecificCommands() {
    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new LoadCommand<DebruijnEnvironment>));

    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new LoadGenomeCommand));

    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new SetMaxVertCommand));
    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new SetFolderCommand));
    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new SetFileNameCommand));

    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new FillPositionCommand));
    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new ClearPositionCommand));

    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new DrawVertexCommand));
    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new DrawEdgeCommand));
    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new DrawPositionCommand));
    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new DrawPartOfGenomeCommand));
    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new DrawContigCommand));
    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new ShowPositionCommand));

    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new PrintPathsCommand));
    AddCommand(shared_ptr<Command<DebruijnEnvironment> >(new PrintContigsStatsCommand));
  }

 public:
  DebruijnOnlineVisualizer() : OnlineVisualizer<DebruijnEnvironment>() {
  }
};

}
