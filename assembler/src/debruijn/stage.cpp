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

void AssemblyStage::load(debruijn_graph::conj_graph_pack& gp) {
    std::string p = path::append_path(cfg::get().output_saves, id_);
    INFO("Loading current state from " << p);

    ScanAll(p, gp, false);
    debruijn_graph::load_lib_data(p);
}

void AssemblyStage::save(const debruijn_graph::conj_graph_pack& gp) const {
    std::string p = path::append_path(cfg::get().output_saves, id_);
    INFO("Saving current state to " << p);

    PrintAll(p, gp);
    debruijn_graph::write_lib_data(p);
}

class StageIdComparator {
  public:
    StageIdComparator(const char* id)
            : id_(id) {}

    bool operator()(const std::unique_ptr<AssemblyStage> &stage) const {
        return 0 == strcmp(id_, stage->id());
    }

  private:
    const char* id_;
};



void StageManager::run(debruijn_graph::conj_graph_pack& g,
                       const char* start_from) {
    auto start_stage = stages_.begin();
    if (start_from) {
        start_stage = std::find_if(stages_.begin(), stages_.end(), StageIdComparator(start_from));
        if (start_stage == stages_.end())
            ERROR("Invalid start stage specified: " << start_from);
        if (start_stage != stages_.begin())
            (*std::prev(start_stage))->load(g);
    }

    for (auto et = stages_.end(); start_stage != et; ++start_stage) {
        AssemblyStage *stage = start_stage->get();

        INFO("STAGE == " << stage->name());
        stage->run(g);
        if (cfg::get().developer_mode && cfg::get().make_saves)
            stage->save(g);
    }

    // For informing spades.py about estimated params
    debruijn_graph::write_lib_data(cfg::get().output_dir);
}

}
