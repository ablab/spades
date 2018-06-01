//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "pipeline/stage.hpp"
#include "io/binary/graphio.hpp"

#include "utils/logger/log_writers.hpp"

#include <algorithm>
#include <cstring>

namespace spades {

const char * const SavesPolicy::CHECKPOINT_FILE = "checkpoint.dat";

void AssemblyStage::load(debruijn_graph::conj_graph_pack& gp,
                         const std::string &load_from,
                         const char* prefix) {
    if (!prefix) prefix = id_;
    std::string p = fs::append_path(load_from, prefix);
    INFO("Loading current state from " << p);

    debruijn_graph::graphio::ScanAll(fs::append_path(p, "graph_pack"), gp);
    debruijn_graph::config::load_lib_data(p);
}

void AssemblyStage::save(const debruijn_graph::conj_graph_pack& gp,
                         const std::string &save_to,
                         const char* prefix) const {
    if (!prefix) prefix = id_;
    std::string p = fs::append_path(save_to, prefix);
    INFO("Saving current state to " << p);
    fs::make_dir(p);

    debruijn_graph::graphio::PrintAll(fs::append_path(p, "graph_pack"), gp);
    debruijn_graph::config::write_lib_data(p);
}

class StageIdComparator {
  public:
    StageIdComparator(const char* id)
            : id_(id) {
        const char* pos = strstr(id, ":");
        len_ = (pos != NULL ? pos - id : strlen(id));
    }

    bool operator()(const std::unique_ptr<AssemblyStage> &stage) const {
        const char* sid = stage->id();
        return (0 == strncmp(id_, sid, len_) && sid[len_] == 0);
    }

  private:
    const char* id_;
    size_t len_;
};

class PhaseIdComparator {
  public:
    PhaseIdComparator(const char* id) {
        const char* pos = strstr(id, ":");
        VERIFY(pos != NULL);
        id_ = pos + 1;
    }

    bool operator()(const std::unique_ptr<CompositeStageBase::PhaseBase> &phase) const {
        return 0 == strcmp(id_, phase->id());
    }

  private:
    const char* id_;
};

void CompositeStageBase::run(debruijn_graph::conj_graph_pack& gp,
                             const char* started_from) {
    // The logic here is as follows. By this time StageManager already called
    // load() function of the Stage itself. Therefore we only need to do
    // storage-related things (if any) and therefore just call the init()
    // function. Phases are supposed only to load the differences.
    VERIFY(parent_);
    init(gp, started_from);
    auto start_phase = phases_.begin();
    if (started_from &&
        strstr(started_from, ":") &&
        started_from == strstr(started_from, id())) {
        start_phase = std::find_if(phases_.begin(), phases_.end(), PhaseIdComparator(started_from));
        if (start_phase == phases_.end()) {
            ERROR("Invalid start stage / phase combination specified: " << started_from);
            exit(-1);
        }
        if (start_phase != phases_.begin()) {
            PhaseBase * prev_phase = std::prev(start_phase)->get();
            std::string composite_id(id());
            composite_id += ":";
            composite_id += prev_phase->id();
            prev_phase->load(gp, parent_->saves_policy().SavesPath(), composite_id.c_str());

        }
    }

    for (auto et = phases_.end(); start_phase != et; ++start_phase) {
        PhaseBase *phase = start_phase->get();

        INFO("PROCEDURE == " << phase->name());
        phase->run(gp, started_from);

        if (parent_->saves_policy().EnabledCheckpoints() != SavesPolicy::Checkpoints::None) {
            std::string composite_id(id());
            composite_id += ":";
            composite_id += phase->id();

            phase->save(gp, parent_->saves_policy().SavesPath(), composite_id.c_str());
            //TODO: erase the previous saves when SavesPolicy::Last
        }
    }

    fini(gp);
}

void StageManager::run(debruijn_graph::conj_graph_pack& g,
                       const char* start_from) {
    auto start_stage = stages_.begin();
    if (start_from) {
        //TODO: refactor
        if (strcmp(start_from, "last") == 0) {
            auto last_saves = saves_policy_.GetLastCheckpoint();
            if (!last_saves.empty()) {
                auto last_stage = std::find_if(stages_.begin(), stages_.end(), StageIdComparator(last_saves.c_str()));
                if (last_stage == stages_.end()) {
                    WARN("Nothing to continue");
                    return;
                }
                start_stage = ++last_stage;
            } else {
                WARN("No saved checkpoint");
            }
        } else {
            start_stage = std::find_if(stages_.begin(), stages_.end(), StageIdComparator(start_from));
            if (start_stage == stages_.end()) {
                ERROR("Invalid start stage specified: " << start_from);
                exit(-1);
            }
        }
        if (start_stage != stages_.begin())
            (*std::prev(start_stage))->load(g, saves_policy_.SavesPath());
    }

    for (; start_stage != stages_.end(); ++start_stage) {
        AssemblyStage *stage = start_stage->get();

        INFO("STAGE == " << stage->name());
        stage->run(g, start_from);
        if (saves_policy_.EnabledCheckpoints() != SavesPolicy::Checkpoints::None) {
            auto prev_saves = saves_policy_.GetLastCheckpoint();
            stage->save(g, saves_policy_.SavesPath());
            saves_policy_.UpdateCheckpoint(stage->id());
            if (!prev_saves.empty() && saves_policy_.EnabledCheckpoints() == SavesPolicy::Checkpoints::Last) {
                fs::remove_dir(fs::append_path(saves_policy_.SavesPath(), prev_saves));
            }
        }
    }
}

}
