//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "mpi_stage.hpp"
#include "partask_mpi.hpp"

#include "pipeline/stage.hpp"

#include "io/binary/graph_pack.hpp"
#include "io/dataset_support/read_converter.hpp"
#include "utils/logger/log_writers.hpp"

#include <algorithm>
#include <vector>
#include <cstring>
#include <mpi.h>


namespace {
class PhaseIdComparator {
  public:
    PhaseIdComparator(const char* id) {
        const char* pos = strstr(id, ":");
        VERIFY(pos != NULL);
        id_ = pos + 1;
    }

    bool operator()(const std::unique_ptr<spades_mpi::MPICompositeStageBase::PhaseBase> &phase) const {
        return 0 == strcmp(id_, phase->id());
    }

  private:
    const char* id_;
};
}

namespace spades_mpi {

void MPICompositeStageBase::run(graph_pack::GraphPack& gp,
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


    // Whether the previous phase was parallel. If this is the first phase then
    // assume that the previous was parallel for the sake of simplicity of the
    // implementation.
    bool pparallel = true;
    for (auto et = phases_.end(); start_phase != et; ++start_phase) {
        PhaseBase *phase = start_phase->get();
        bool cparallel = phase->distributed();

        if (cparallel) {
            if (!pparallel) {
                partask::critical_ordered([this] {
                    if (worker()) {
                        io::ConvertIfNeeded(cfg::get_writable().ds.reads, cfg::get().max_threads);
                    }
                });
                INFO("Syncing world for MPI parallel section");
                const size_t deadbeef = 0xDEADBEEF;
                if (master()) {
                    partask::OutputMPIStreamBcast s(0);
                    io::binary::FullPackIO().BinWrite(s, gp);
                    io::binary::BinWrite(s, deadbeef);
                    debruijn_graph::config::write_lib_data(s);
                    io::binary::BinWrite(s, deadbeef);
                } else {
                    partask::InputMPIStreamBcast s(0);
                    io::binary::FullPackIO().BinRead(s, gp);
                    size_t db;
                    io::binary::BinRead(s, db);
                    VERIFY(db == deadbeef);
                    debruijn_graph::config::load_lib_data(s);
                    io::binary::BinRead(s, db);
                    VERIFY(db == deadbeef);
                }
                INFO("World synced");
            }
            INFO("MPI PROCEDURE == " << phase->name() << (master() ? " (master)" : " (worker)"));
            phase->run(gp, started_from);

            // Do saves only on master node
            if (parent_->saves_policy().EnabledCheckpoints(id()) && master()) {
                std::string composite_id(id());
                composite_id += ":";
                composite_id += phase->id();

                phase->save(gp, parent_->saves_policy().SavesPath(), composite_id.c_str());
            }
        } else {
            if (master()) {
                INFO("PROCEDURE == " << phase->name());
                phase->run(gp, started_from);
                if (parent_->saves_policy().EnabledCheckpoints(id())) {
                    std::string composite_id(id());
                    composite_id += ":";
                    composite_id += phase->id();

                    phase->save(gp, parent_->saves_policy().SavesPath(), composite_id.c_str());
                }
            } else {
                INFO("PROCEDURE == " << phase->name() << " (skipped on worker)");
            }
        }

        pparallel = cparallel;
    }

    fini(gp);
}

MPIStageManager::MPIStageManager(spades::SavesPolicy policy)
        : StageManager(policy), world_size_(1), rank_(0), first_(false) {
    int initialized = 0;
    MPI_Initialized(&initialized);
    VERIFY(initialized);
    if (!initialized) {
        int provided;
        MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, &provided);
        if (provided < MPI_THREAD_FUNNELED) {
            FATAL_ERROR("Used MPI implementation failed to provide MPI_THREAD_FUNNELED thread support level");
        }
        first_ = true;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &world_size_);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

    INFO("MPI communications established, world size: " << world_size_ << ", current rank: " << rank_ << (master() ? " (master)" : " (worker)"));
}

MPIStageManager::~MPIStageManager() {
    INFO("MPI communications stopped" << (master() ? " (master)" : " (worker)"));
    if (first_)
        MPI_Finalize();
}

void MPIStageManager::run(graph_pack::GraphPack& g,
                          const char* start_from) {
    auto start_stage = prepare_run(g, start_from);
    const auto& saves_policy = this->saves_policy();

    // Whether the previous stage was parallel. If this is the first stage then
    // assume that the previous was parallel for the sake of simplicity of the
    // implementation.
    bool pparallel = true;

    for (auto current_stage = stages().begin(); current_stage != stages().end(); ++current_stage) {
        spades::AssemblyStage *stage = current_stage->get();
        if (current_stage < start_stage && !stage->run_on_load()) {
            continue;
        }

        bool cparallel = stage->distributed();

        if (cparallel) {
            if (!pparallel) {
                partask::critical_ordered([this] {
                    if (worker()) {
                        io::ConvertIfNeeded(cfg::get_writable().ds.reads, cfg::get().max_threads);
                    }
                });
                INFO("Syncing world for MPI parallel section");
                const size_t deadbeef = 0xDEADBEEF;
                if (master()) {
                    partask::OutputMPIStreamBcast s(0);
                    io::binary::FullPackIO().BinWrite(s, g);
                    io::binary::BinWrite(s, deadbeef);
                    debruijn_graph::config::write_lib_data(s);
                    io::binary::BinWrite(s, deadbeef);
                } else {
                    partask::InputMPIStreamBcast s(0);
                    io::binary::FullPackIO().BinRead(s, g);
                    size_t db;
                    io::binary::BinRead(s, db);
                    VERIFY_MSG(db == deadbeef, "Values " << db << " " << deadbeef);
                    debruijn_graph::config::load_lib_data(s);
                    io::binary::BinRead(s, db);
                    VERIFY(db == deadbeef);
                }
                INFO("World synced");
            }
            INFO("MPI STAGE == " << stage->name() << (master() ? " (master)" : " (worker)"));
            stage->prepare(g, start_from);
            stage->run(g, start_from);

            // Do saves only on master node
            if (saves_policy.EnabledCheckpoints(stage->id()) && master())
                stage->save(g, saves_policy.SavesPath());
        } else {
            if (master()) {
                INFO("STAGE == " << stage->name());
                stage->prepare(g, start_from);
                stage->run(g, start_from);
                if (saves_policy.EnabledCheckpoints(stage->id())) {
                    auto prev_saves = saves_policy.GetLastCheckpoint();
                    stage->save(g, saves_policy.SavesPath());
                    saves_policy.UpdateCheckpoint(stage->id());
                    if (!prev_saves.empty() && saves_policy.RemovePreviousCheckpoint()) {
                        remove_all(saves_policy.SavesPath() / prev_saves);
                    }
                }
            } else {
                INFO("STAGE == " << stage->name() << " (skipped on worker)");
            }
        }

        if (cparallel || !stage->constant()) {
            pparallel = cparallel;
        }
    }
}

bool MPIAssemblyStage::master() const { return static_cast<const MPIStageManager*>(parent_)->master(); }
bool MPIAssemblyStage::worker() const { return static_cast<const MPIStageManager*>(parent_)->worker(); }

}  // namespace spades
