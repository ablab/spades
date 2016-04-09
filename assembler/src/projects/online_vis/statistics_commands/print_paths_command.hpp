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

class PrintPathsCommand : public LocalCommand<DebruijnEnvironment> {
  typedef vector<EdgeId> Path;

 protected:
  size_t MinArgNumber() const {
    return 2;
  }

  bool CheckCorrectness(const vector<string>& args) const {
    return CheckEnoughArguments(args);
  }

 public:
  string Usage() const {
    string answer;
    answer = answer + "Command `paths` \n" +
      "Usage:\n" +
      "> paths <vertex_from> <vertex_to> [<max_length>] \n" +
      " This command prints all paths between two given vertices, that do not exceed `max_length` parameter.\n" +
      " You should specify two integers (id of vertices), between which you want to find paths." +
      " Optionally you can provide `max_length` integer, \n" +
      " so that tool does not consider paths longer than `max_length`. It is equal to 100000 by default.";
    return answer;
  }

  PrintPathsCommand() : LocalCommand<DebruijnEnvironment>("print_paths")
  {
  }

  void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
    const vector<string>& args = arg_list.GetAllArguments();
    if (!CheckCorrectness(args))
      return;

    TRACE("Executing `paths` command");
//    bool first_edge = false;
//    bool second_edge = false;
    size_t from = GetInt(args[1]);
    size_t to = GetInt(args[2]);
    size_t max_length = 100000;
    if (args.size() > 3)
      max_length = GetInt(args[3]);
    if (arg_list["edge"]=="true") {
      // not so good, maybe do something
//      first_edge = second_edge = true;
      TRACE("Looking at edges");
      if (!CheckEdgeExists(curr_env.finder(), from) || !CheckEdgeExists(curr_env.finder(), to))
        return;
      from = curr_env.graph().int_id(curr_env.graph().EdgeEnd(curr_env.finder().ReturnEdgeId(from)));
      to = curr_env.graph().int_id(curr_env.graph().EdgeStart(curr_env.finder().ReturnEdgeId(to)));
    }

    if (!CheckVertexExists(curr_env.finder(), from) || !CheckVertexExists(curr_env.finder(), to))
      return;
  
    const Graph& graph = curr_env.graph();

    TRACE("Looking for the paths");
    PathStorageCallback<Graph> callback(curr_env.graph());
    ProcessPaths(graph, 0, max_length,
        curr_env.finder().ReturnVertexId(from),
        curr_env.finder().ReturnVertexId(to), callback);
    const vector<Path>& paths = callback.paths();

    cout << paths.size() << " path(s) have been found : " << endl;
    for (size_t i = 0; i < paths.size(); ++i) {
      cout << (i + 1) << "-th path (" << PathLength(graph, paths[i]) << ") ::: ";
      for (size_t j = 0; j < paths[i].size(); ++j) {
        cout << paths[i][j].int_id()
             << "(" << graph.length(paths[i][j]) << ") ";
      }
      cout << endl;
    }

  }

 private:
  size_t PathLength(const Graph& g, const vector<EdgeId>& path) const {
    size_t res = 0;
    for (size_t i = 0; i < path.size(); ++i)
      res += g.length(path[i]);
    return res;
  }
  
  DECL_LOGGER("PrintPathsCommand");
};

}
