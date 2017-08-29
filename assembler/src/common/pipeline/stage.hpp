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
            : name_(name), id_(id), parent_(nullptr) { }

    virtual ~AssemblyStage() = default;
    AssemblyStage(const AssemblyStage &) = delete;
    AssemblyStage &operator=(const AssemblyStage &) = delete;

    const char *name() const { return name_; }
    const char *id() const { return id_; }

    virtual void load(debruijn_graph::conj_graph_pack &, const std::string &load_from, const char *prefix = nullptr);
    virtual void save(const debruijn_graph::conj_graph_pack &, const std::string &save_to,
                      const char *prefix = nullptr) const;
    virtual void run(debruijn_graph::conj_graph_pack &, const char *started_from = nullptr) = 0;

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
                : AssemblyStage(name, id), parent_stage_(nullptr) { }

    protected:
        CompositeStageBase *parent_stage_;

        friend class CompositeStageBase;
    };

    CompositeStageBase(const char *name, const char *id)
            : AssemblyStage(name, id) { }

    CompositeStageBase &add(PhaseBase *phase) {
        phases_.push_back(std::unique_ptr<PhaseBase>(phase));
        phase->parent_stage_ = this;

        return *this;
    }

    template<typename Phase, typename ... Args>
    CompositeStageBase &add(Args&&... args) {
        phases_.push_back(std::unique_ptr<Phase>(new Phase(std::forward<Args>(args)...)));
        phases_.back()->parent_stage_ = this;

        return *this;
    }

    virtual void init(debruijn_graph::conj_graph_pack &, const char * = nullptr) = 0;
    virtual void fini(debruijn_graph::conj_graph_pack &) = 0;
    void run(debruijn_graph::conj_graph_pack &gp, const char * = nullptr);

private:
    std::vector<std::unique_ptr<PhaseBase> > phases_;
};

template<class Storage>
class CompositeStageWithStorage : public CompositeStageBase {
public:
    class Phase : public PhaseBase {
    public:
        Phase(const char *name, const char *id)
                : PhaseBase(name, id) { }

        CompositeStageWithStorage<Storage> *parent() { return static_cast<CompositeStageWithStorage<Storage> *>(parent_stage_); }
        const CompositeStageWithStorage<Storage> *parent() const { return static_cast<const CompositeStageWithStorage<Storage> *>(parent_stage_); }

        Storage &storage() { return parent()->storage(); }
        const Storage &storage() const { return parent()->storage(); }
    };

    CompositeStageWithStorage(const char *name, const char *id)
            : CompositeStageBase(name, id) { }

    void init(debruijn_graph::conj_graph_pack &, const char * = nullptr) override {};
    void fini(debruijn_graph::conj_graph_pack &) override {};

    virtual Storage &storage() = 0;
    virtual const Storage &storage() const = 0;
};

// FIXME: Make storage a policy
template<class Storage>
class CompositeStage : public CompositeStageWithStorage<Storage> {
public:
    CompositeStage(const char *name, const char *id)
            : CompositeStageWithStorage<Storage>(name, id) { }

    Storage &storage() override { return storage_; }
    const Storage &storage() const override { return storage_; }

private:
    Storage storage_;
};

template<class Storage>
class CompositeStageDeferred : public CompositeStageWithStorage<Storage> {
public:
    CompositeStageDeferred(const char *name, const char *id)
            : CompositeStageWithStorage<Storage>(name, id) { }

    Storage &storage() override { return *storage_; }
    const Storage &storage() const override { return *storage_; }

protected:
    bool has_storage() const { return (bool)storage_; }

    template <typename...Args> void init_storage(Args&&... args) {
        storage_.reset(new Storage(std::forward<Args>(args)...));
    }
    void reset_storage() {
        storage_.reset();
    }

private:
    // std::optional would be better, but it requires complete Storage type at
    // this point.
    std::unique_ptr<Storage> storage_;
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

    template<typename Stage, typename ... Args>
    StageManager &add(Args&&... args) {
        stages_.push_back(std::unique_ptr<Stage>(new Stage(std::forward<Args>(args)...)));
        stages_.back()->parent_ = this;

        return *this;
    }

    void run(debruijn_graph::conj_graph_pack &g,
             const char *start_from = nullptr);

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
