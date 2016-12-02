//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "visualization/position_filler.hpp"
#include "pipeline/graph_pack.hpp"
#include "visualization/visualization_utils.hpp"
#include "standard_vis.hpp"

namespace online_visualization {

typedef debruijn_graph::BasicSequenceMapper<debruijn_graph::Graph, Index> MapperClass;
typedef visualization::position_filler::PosFiller<Graph> FillerClass;
typedef debruijn_graph::KmerMapper<Graph> KmerMapperClass;
typedef omnigraph::GraphElementFinder<Graph> ElementFinder;
typedef shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> ColoringClass;

class Environment : private boost::noncopyable {
 protected:
  const string name_;
  const string path_;

 public:
  Environment(const string &name, const string &path)
      : name_(name),
        path_(path) {
  }

  virtual ~Environment() {
  }

  inline string name() const {
    return name_;
  }

  inline string path() const {
    return path_;
  }

  virtual string str() const {
    stringstream ss;
    ss << name_ + " " + path_;
    return ss.str();
  }

  virtual inline bool IsCorrect() const {
    // make here some checks! for path etc
    return true;
  }

};

}
