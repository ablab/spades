//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "graph_pack.hpp"
#include "omni/visualization_utils.hpp"
#include "standard_vis.hpp"
#include "debruijn_stats.hpp"

namespace online_visualization {

typedef debruijn_graph::NewExtendedSequenceMapper<Graph> MapperClass;
typedef debruijn_graph::PosFiller<Graph, MapperClass> FillerClass;
typedef debruijn_graph::KmerMapper<Graph> KmerMapperClass;
typedef map<EdgeId, string> ColoringClass;

class Environment {
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
