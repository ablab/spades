//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

#include <vector>
#include <memory>
#include <type_traits>

namespace spades {

class MPIAssemblyStage;

class MPIAssemblyStage : public AssemblyStage {
public:
    using AssemblyStage::AssemblyStage;

    bool master() const;
    bool worker() const;
    bool distributed() const override { return true; }
};

class MPICompositeStageBase : public MPIAssemblyStage {
public:
    class PhaseBase : public MPIAssemblyStage {
    public:
        PhaseBase(const char *name, const char *id)
                : MPIAssemblyStage(name, id), parent_stage_(nullptr) { }

        bool distributed() const override { return false; }
        bool master() const { return parent_stage_->master(); }
        bool worker() const { return parent_stage_->worker(); }
    protected:
        MPICompositeStageBase *parent_stage_;

        friend class MPICompositeStageBase;
    };

    MPICompositeStageBase(const char *name, const char *id)
            : MPIAssemblyStage(name, id) { }

    MPICompositeStageBase &add(PhaseBase *phase) {
        phases_.push_back(std::unique_ptr<PhaseBase>(phase));
        phase->parent_stage_ = this;

        return *this;
    }

    template<typename Phase, typename ... Args>
    MPICompositeStageBase &add(Args&&... args) {
        phases_.push_back(std::unique_ptr<Phase>(new Phase(std::forward<Args>(args)...)));
        phases_.back()->parent_stage_ = this;

        return *this;
    }

    const std::vector<std::unique_ptr<PhaseBase> >& phases() const {
        return phases_;
    }

    virtual void init(graph_pack::GraphPack &, const char * = nullptr) = 0;
    virtual void fini(graph_pack::GraphPack &) = 0;
    void run(graph_pack::GraphPack &gp, const char * = nullptr);

private:
    std::vector<std::unique_ptr<PhaseBase> > phases_;
};

template<class Storage>
class MPICompositeStageWithStorage : public MPICompositeStageBase {
public:
    class Phase : public PhaseBase {
    public:
        Phase(const char *name, const char *id)
                : PhaseBase(name, id) { }

        MPICompositeStageWithStorage<Storage> *parent() { return static_cast<MPICompositeStageWithStorage<Storage> *>(parent_stage_); }
        const MPICompositeStageWithStorage<Storage> *parent() const { return static_cast<const MPICompositeStageWithStorage<Storage> *>(parent_stage_); }

        Storage &storage() { return parent()->storage(); }
        const Storage &storage() const { return parent()->storage(); }
    };

    MPICompositeStageWithStorage(const char *name, const char *id)
            : MPICompositeStageBase(name, id) { }

    void init(graph_pack::GraphPack &, const char * = nullptr) override {};
    void fini(graph_pack::GraphPack &) override {};

    virtual Storage &storage() = 0;
    virtual const Storage &storage() const = 0;
};

// FIXME: Make storage a policy
template<class Storage>
class MPICompositeStage : public MPICompositeStageWithStorage<Storage> {
public:
    MPICompositeStage(const char *name, const char *id)
            : MPICompositeStageWithStorage<Storage>(name, id) { }

    Storage &storage() override { return storage_; }
    const Storage &storage() const override { return storage_; }

private:
    Storage storage_;
};

template<class Storage>
class MPICompositeStageDeferred : public MPICompositeStageWithStorage<Storage> {
public:
    MPICompositeStageDeferred(const char *name, const char *id)
            : MPICompositeStageWithStorage<Storage>(name, id) { }

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

class MPIStageManager : public StageManager {
public:
    MPIStageManager(SavesPolicy policy = SavesPolicy());
    ~MPIStageManager();

    void run(graph_pack::GraphPack &g,
             const char *start_from = nullptr) override;

    bool master() const { return rank_ == 0; }
    bool worker() const { return rank_ != 0; }

private:
    int world_size_;
    int rank_;
    bool first_;

    DECL_LOGGER("MPIStageManager");
};

}  // namespace spades
