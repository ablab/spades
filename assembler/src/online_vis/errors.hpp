//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "vis_utils.hpp"

namespace online_visualization {

  void FireGenericError(const string& msg) {
    cout << msg << endl;
    cout << "Please try again" << endl;
  }

  void FireEdgeDoesNotExist(size_t edge_id) {
    cout << "Ignoring the request. The edge " << edge_id << " does not exist" << endl;
    cout << "Please try again" << endl;
  }

  void FireVertexDoesNotExist(size_t vertex_id) {
    cout << "Ignoring the request. The vertex " << vertex_id << " does not exist" << endl;
    cout << "Please try again" << endl;
  }

  void FireNoCorrespondingGraphLocation(string location) {
    cout << "No corresponding graph location " << location << endl;
  }

  void FireNotEnoughArguments() {
    cout << "Not enough arguments" << endl;
    cout << "Please try again" << endl;
  }

  void FireFileDoesNotExist(const string& file) {
    cout << "File " << file << " does not exist." << endl;
    cout << "Please try again" << endl;
  }

  void FireBadArgument(const string& arg) {
    cout << "Bad word specifier: `" << arg << "'" << endl;
    cout << "Please try again" << endl;
  }

  void FireNumberOutOfBounds(int num_of_command) {
    cout << "The command number parameter " << num_of_command
         << " must be positive and not exceed the size of history" << endl;
    cout << "Please try again" << endl;
  }

  bool CheckFileExists(const string& file) {
    if (!path::is_regular_file(file)) {
      FireFileDoesNotExist(file);
      return false;
    }
    return true;
  }

  bool CheckPositionBounds(size_t position, size_t total_size, size_t K) {
    bool result = (position + K + 1) <= total_size;
    if (!result) {
      cout << "Ignoring the request. Position is out of range : required position is "
        << position << " while length of the sequence is "
        << total_size << endl;
      cout << "Please try again" << endl;
    }
    return result;
  }

  bool CheckIsNumber(const string& str) {
    if (!IsNumber(str)) {
      cout << "The argument `" << str << "' is not a number" << endl;
      cout << "Please try again" << endl;
      return false;
    }
    return true;
  }

  bool CheckEnvIsCorrect(string path, size_t K) {
    if (!CheckFileExists(path + ".grp"))
      return false;
    if (!CheckFileExists(path + ".sqn"))
      return false;

    if (!(K >= runtime_k::MIN_K && cfg::get().K < runtime_k::MAX_K)) {
      LOG("K " << K << " is out of bounds");
      return false;
    }
    if (K % 2 == 0) {
      LOG("K must be odd");
      return false;
    }

    return true;
  }

  bool CheckVertexExists(const GraphElementFinder<Graph>& finder, size_t vertex_id) {
    VertexId vertex = finder.ReturnVertexId(vertex_id);
    if (vertex == VertexId(NULL)) {
      FireVertexDoesNotExist(vertex_id);
      return false;
    }
    else {
      return true;
    }
  }

  bool CheckEdgeExists(const GraphElementFinder<Graph>& finder, size_t edge_id) {
    EdgeId edge = finder.ReturnEdgeId(edge_id);
    if (edge == EdgeId(NULL)) {
      FireEdgeDoesNotExist(edge_id);
      return false;
    }
    else {
      return true;
    }
  }

}
