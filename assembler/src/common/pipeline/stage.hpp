//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __STAGE_HPP__
#define __STAGE_HPP__

#include "pipeline/graph_pack.hpp"

#include <vector>
#include <memory>

namespace spades {

class StageManager;

class AssemblyStage {
public:
    AssemblyStage(const char *name, const char *id)
            : name_(name), id_(id), parent_(NULL) { }

    virtual ~AssemblyStage() { }

    AssemblyStage(const AssemblyStage &) = delete;

    AssemblyStage &operator=(const AssemblyStage &) = delete;

    const char *name() const { return name_; }

    const char *id() const { return id_; }

    virtual void load(debruijn_graph::conj_graph_pack &, const std::string &load_from, const char *prefix = NULL);

    virtual void save(const debruijn_graph::conj_graph_pack &, const std::string &save_to,
                      const char *prefix = NULL) const;

    virtual void run(debruijn_graph::conj_graph_pack &, const char *started_from = NULL) = 0;

private:
    const char *name_;
    const char *id_;

protected:
    const StageManager *parent_;

    friend class StageManager;
};

class CompositeStageBase : public AssemblyStage {
public:
    class PhaseBase : public AssemblyStage {
    public:
        PhaseBase(const char *name, const char *id)
                : AssemblyStage(name, id), parent_stage_(NULL) { }

    protected:
        CompositeStageBase *parent_stage_;

        friend class CompositeStageBase;
    };

    CompositeStageBase(const char *name, const char *id)
            : AssemblyStage(name, id) { }

    CompositeStageBase *add(PhaseBase *phase) {
        phases_.push_back(std::unique_ptr<PhaseBase>(phase));
        phase->parent_stage_ = this;

        return this;
    }

    CompositeStageBase *add(std::initializer_list<PhaseBase *> phases) {
        for (auto it = phases.begin(), et = phases.end(); it != et; ++it)
            add(*it);

        return this;
    }

    void run(debruijn_graph::conj_graph_pack &gp, const char * = NULL);

private:
    std::vector<std::unique_ptr<PhaseBase> > phases_;
};

template<class Storage>
class CompositeStage : public CompositeStageBase {
public:
    class Phase : public PhaseBase {
    public:
        Phase(const char *name, const char *id)
                : PhaseBase(name, id) { }

        CompositeStage<Storage> *parent() { return static_cast<CompositeStage<Storage> *>(parent_stage_); }

        const CompositeStage<Storage> *parent() const { return static_cast<const CompositeStage<Storage> *>(parent_stage_); }

        Storage &storage() { return parent()->storage(); }

        const Storage &storage() const { return parent()->storage(); }
    };

    CompositeStage(const char *name, const char *id)
            : CompositeStageBase(name, id) { }

    Storage &storage() { return storage_; }

    const Storage &storage() const { return storage_; }

private:
    Storage storage_;
};

class StageManager {

public:
    struct SavesPolicy {
        bool make_saves_;
        std::string load_from_;
        std::string save_to_;

        SavesPolicy()
                : make_saves_(false), load_from_(""), save_to_("") { }

        SavesPolicy(bool make_saves, const std::string &load_from, const std::string &save_to)
                : make_saves_(make_saves), load_from_(load_from), save_to_(save_to) { }
    };

    StageManager(SavesPolicy policy = SavesPolicy())
            : saves_policy_(policy) { }

    StageManager &add(AssemblyStage *stage) {
        stages_.push_back(std::unique_ptr<AssemblyStage>(stage));
        stages_.back()->parent_ = this;

        return *this;
    }

    StageManager &add(std::initializer_list<AssemblyStage *> stages) {
        for (auto it = stages.begin(), et = stages.end(); it != et; ++it)
            add(*it);

        return *this;
    }

    void run(debruijn_graph::conj_graph_pack &g,
             const char *start_from = NULL);

    const SavesPolicy &saves_policy() const {
        return saves_policy_;
    }

private:
    std::vector<std::unique_ptr<AssemblyStage> > stages_;
    SavesPolicy saves_policy_;

    DECL_LOGGER("StageManager");
};


};

#endif // __STAGE_HPP__
