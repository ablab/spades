//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef __STAGE_HPP__
#define __STAGE_HPP__

#include "graph_pack.hpp"

#include <vector>
#include <memory>

namespace spades {

class AssemblyStage {
  public:
    AssemblyStage(const char *name, const char *id)
            : name_(name), id_(id) {}

    AssemblyStage(const AssemblyStage&) = delete;
    AssemblyStage& operator=(const AssemblyStage&) = delete;

    const char *name() const { return name_; }
    const char *id() const { return id_; }

    virtual void load(debruijn_graph::conj_graph_pack&, const char* prefix = NULL);
    virtual void save(const debruijn_graph::conj_graph_pack&, const char* prefix = NULL) const;
    virtual void run(debruijn_graph::conj_graph_pack&, const char* started_from = NULL) = 0;

  private:
    const char *name_;
    const char *id_;
};

class CompositeStageBase : public AssemblyStage {
  public:
    class PhaseBase : public AssemblyStage {
      public:
        PhaseBase(const char *name, const char *id)
                : AssemblyStage(name, id), parent_(NULL) {}
      protected:
        CompositeStageBase *parent_;

        friend class CompositeStageBase;
    };

    CompositeStageBase(const char *name, const char *id)
            : AssemblyStage(name, id) {}

    CompositeStageBase* add(PhaseBase *phase) {
        phases_.push_back(std::unique_ptr<PhaseBase>(phase));
        phase->parent_ = this;

        return this;
    }

    void run(debruijn_graph::conj_graph_pack& gp, const char* = NULL);

  private:
    std::vector<std::unique_ptr<PhaseBase> > phases_;
};

template<class Storage>
class CompositeStage : public CompositeStageBase {
  public:
    class Phase : public PhaseBase {
      public:
        Phase(const char *name, const char *id)
                : PhaseBase(name, id) {}

        CompositeStage<Storage>* parent() { return static_cast<CompositeStage<Storage>*>(parent_); }
        const CompositeStage<Storage>* parent() const { return static_cast<const CompositeStage<Storage>*>(parent_); }

        Storage &storage() { return parent()->storage(); }
        const Storage &storage() const { return parent()->storage(); }
    };

    CompositeStage(const char *name, const char *id)
            : CompositeStageBase(name, id) {}

    Storage &storage() { return storage_; }
    const Storage &storage() const { return storage_; }

  private:
    Storage storage_;
};

class StageManager {
  public:
    void add(AssemblyStage *stage) {
        stages_.push_back(std::unique_ptr<AssemblyStage>(stage));
    }
    void run(debruijn_graph::conj_graph_pack& g,
             const char* start_from = NULL);
  private:
    std::vector<std::unique_ptr<AssemblyStage> > stages_;

    DECL_LOGGER("StageManager");
};


};

#endif // __STAGE_HPP__
