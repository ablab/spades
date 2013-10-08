//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "stage.hpp"
#include "graphio.hpp"

#include "logger/log_writers.hpp"

#include <algorithm>
#include <cstring>

namespace spades {

void AssemblyStage::load(debruijn_graph::conj_graph_pack& gp, const char* prefix) {
    std::string p = path::append_path(cfg::get().load_from, prefix == NULL ? id_ : prefix);
    INFO("Loading current state from " << p);

    ScanAll(p, gp, false);
    debruijn_graph::load_lib_data(p);
}

void AssemblyStage::save(const debruijn_graph::conj_graph_pack& gp, const char* prefix) const {
    std::string p = path::append_path(cfg::get().output_saves, prefix == NULL ? id_ : prefix);
    INFO("Saving current state to " << p);

    PrintAll(p, gp);
    debruijn_graph::write_lib_data(p);
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
            prev_phase->load(gp, composite_id.c_str());
        }
    }

    for (auto et = phases_.end(); start_phase != et; ++start_phase) {
        PhaseBase *phase = start_phase->get();

        INFO("PROCEDURE == " << phase->name());
        phase->run(gp, started_from);

        if (cfg::get().developer_mode && cfg::get().make_saves) {
            std::string composite_id(id());
            composite_id += ":";
            composite_id += phase->id();

            phase->save(gp, composite_id.c_str());
        }

    }
}

void StageManager::run(debruijn_graph::conj_graph_pack& g,
                       const char* start_from) {
    auto start_stage = stages_.begin();
    if (start_from) {
        start_stage = std::find_if(stages_.begin(), stages_.end(), StageIdComparator(start_from));
        if (start_stage == stages_.end()) {
            ERROR("Invalid start stage specified: " << start_from);
            exit(-1);
        }
        if (start_stage != stages_.begin())
            (*std::prev(start_stage))->load(g);
    }

    for (auto et = stages_.end(); start_stage != et; ++start_stage) {
        AssemblyStage *stage = start_stage->get();

        INFO("STAGE == " << stage->name());
        stage->run(g, start_from);
        if (cfg::get().developer_mode && cfg::get().make_saves)
            stage->save(g);
    }

    // For informing spades.py about estimated params
    debruijn_graph::write_lib_data(cfg::get().output_dir);
}

}
