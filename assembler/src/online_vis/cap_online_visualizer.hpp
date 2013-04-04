#pragma once

#include "online_visualizer.hpp"
#include "cap_environment.hpp"
#include "cap_commands.hpp"

namespace online_visualization {

class CapOnlineVisualizer : public OnlineVisualizer<CapEnvironment> {
 protected:
  void AddSpecificCommands() {
    AddCommand(shared_ptr<Command<CapEnvironment> >(new LoadCommand<CapEnvironment>));
    AddCommand(shared_ptr<Command<CapEnvironment> >(new AddGenomeCommand));
    AddCommand(shared_ptr<Command<CapEnvironment> >(new BuildGraphCommand));
    AddCommand(shared_ptr<Command<CapEnvironment> >(new RefineCommand));
    AddCommand(shared_ptr<Command<CapEnvironment> >(new SaveGenomesCommand));
    AddCommand(shared_ptr<Command<CapEnvironment> >(new SaveGraphCommand));
    AddCommand(shared_ptr<Command<CapEnvironment> >(new FindIndelsCommand));
    AddCommand(shared_ptr<Command<CapEnvironment> >(new FindInversionsCommand));
    //AddCommand(shared_ptr<Command<CapEnvironment> >(new LoadGraphCommand));
  }

 public:
  CapOnlineVisualizer() : OnlineVisualizer<CapEnvironment>() {
  }
};

}
