//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef __STAGE_HPP__
#define __STAGE_HPP__

#include "graph_pack.hpp"

#include <vector>

namespace spades {

class AssemblyStage {
  public:
    AssemblyStage(const char *name, const char *id)
            : name_(name), id_(id) {}

    AssemblyStage(const AssemblyStage&) = delete;
    AssemblyStage& operator=(const AssemblyStage&) = delete;

    const char *name() const { return name_; }
    const char *id() const { return id_; }

    virtual void load(debruijn_graph::conj_graph_pack&);
    virtual void save(const debruijn_graph::conj_graph_pack&) const;
    virtual void run(debruijn_graph::conj_graph_pack&) = 0;

  private:
    const char *name_;
    const char *id_;
};

class StageManager {
  public:
    void add(AssemblyStage* stage) {
        stages_.push_back(stage);
    }
    void run(debruijn_graph::conj_graph_pack& g,
             const char* start_from = NULL);
  private:
    std::vector<AssemblyStage*> stages_;

    DECL_LOGGER("StageManager");
};


};

#endif // __STAGE_HPP__
