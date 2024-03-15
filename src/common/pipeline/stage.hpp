//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __STAGE_HPP__
#define __STAGE_HPP__

#include "graph_pack.hpp"

#include "configs/config_struct.hpp"
#include "utils/logger/logger.hpp"

#include <filesystem>
#include <memory>
#include <variant>
#include <vector>

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

    /// @throw std::ios_base::failure if load_from does not contain all required files

    virtual void load(graph_pack::GraphPack &, const std::filesystem::path &load_from, const char *prefix = nullptr);
    virtual void save(const graph_pack::GraphPack &, const std::filesystem::path &save_to,
                      const char *prefix = nullptr) const;
    void prepare(graph_pack::GraphPack &, const char *stage_name, const char *started_from = nullptr);
    virtual void run(graph_pack::GraphPack &, const char *started_from = nullptr) = 0;

private:
    const char *name_;

protected:
    const char *id_;
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

    virtual void init(graph_pack::GraphPack &, const char * = nullptr) = 0;
    virtual void fini(graph_pack::GraphPack &) = 0;
    void run(graph_pack::GraphPack &gp, const char * = nullptr);

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

    void init(graph_pack::GraphPack &, const char * = nullptr) override {};
    void fini(graph_pack::GraphPack &) override {};

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

class SavesPolicy {
public:
    using Checkpoints = debruijn_graph::config::Checkpoints;

    SavesPolicy()
            : checkpoints_(Checkpoints::None), saves_path_("") {
    }

    SavesPolicy(const std::variant<Checkpoints, std::string>& checkpoints,
                const std::filesystem::path &saves_path, const std::filesystem::path &load_path = "")
            : checkpoints_(checkpoints), saves_path_(saves_path) {
        load_path_ = (load_path == "" ? saves_path_ : load_path);
    }

    bool EnabledAnyCheckpoint() const {
        return std::holds_alternative<std::string>(checkpoints_) or std::get<Checkpoints>(checkpoints_) != SavesPolicy::Checkpoints::None;
    }

    bool EnabledCheckpoints(const std::string& stage_id) const {
        if (std::holds_alternative<std::string>(checkpoints_))
            return stage_id.compare(std::get<std::string>(checkpoints_)) == 0;
        return std::get<Checkpoints>(checkpoints_) != SavesPolicy::Checkpoints::None;
    }

    bool RemovePreviousCheckpoint() const {
        return std::holds_alternative<Checkpoints>(checkpoints_) and std::get<Checkpoints>(checkpoints_) == SavesPolicy::Checkpoints::Last;
    }

    const std::filesystem::path & SavesPath() const { return saves_path_; }
    const std::filesystem::path & LoadPath() const { return load_path_; }

    std::string GetLastCheckpoint() const {
        std::string res;
        std::ifstream ifs(saves_path_ / CHECKPOINT_FILE);
        if (ifs.is_open())
            ifs >> res;
        return res;
    }

    void UpdateCheckpoint(const char *name) const {
        std::ofstream(saves_path_ / CHECKPOINT_FILE) << name;
    }

private:
    static constexpr const char *CHECKPOINT_FILE = "checkpoint.dat";

    std::variant<Checkpoints, std::string> checkpoints_;
    std::filesystem::path saves_path_;
    std::filesystem::path load_path_;
};

class StageManager {
public:
    StageManager(SavesPolicy policy = SavesPolicy())
            : saves_policy_(std::move(policy)) { }
    virtual ~StageManager() = default;

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

    void run(graph_pack::GraphPack &g,
             const char *start_from = nullptr);

    const SavesPolicy &saves_policy() const {
        return saves_policy_;
    }

private:
    using Stages = std::vector<std::unique_ptr<AssemblyStage> >;

    Stages stages_;
    SavesPolicy saves_policy_;

    DECL_LOGGER("StageManager");
};


};

#endif // __STAGE_HPP__
