//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../errors.hpp"

namespace online_visualization {

class PrintEdgeCommand : public LocalCommand<DebruijnEnvironment> {

 protected:
  size_t MinArgNumber() const {
    return 1;
  }

  bool CheckCorrectness(const vector<string>& args) const {
    return CheckEnoughArguments(args);
  }

 public:
  string Usage() const {
    string answer;
    answer = answer + "Command `paths` \n" +
      "Usage:\n" +
      "> print_edge <edge_id> \n" +
      " This command prints edge coverage and sequence.";
    return answer;
  }

  PrintEdgeCommand() : LocalCommand<DebruijnEnvironment>("print_edge")
  {
  }

  void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
    const vector<string>& args = arg_list.GetAllArguments();
    if (!CheckCorrectness(args))
      return;
    TRACE("Executing `print_edge` command");
    size_t edgeID = GetInt(args[1]);
    if (!CheckEdgeExists(curr_env.finder(), edgeID))
        return;
    EdgeId edge = curr_env.finder().ReturnEdgeId(edgeID);
    cout << curr_env.graph().str(edge) << endl;

    cout << curr_env.graph().EdgeNucls(edge) << endl;


  }

  DECL_LOGGER("PrintEdgeCommand");
};

}
