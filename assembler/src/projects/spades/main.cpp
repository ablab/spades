//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "pipeline/config_struct.hpp"

#include "utils/logger/log_writers.hpp"
#include "utils/memory_limit.hpp"
#include "utils/segfault_handler.hpp"
#include "utils/filesystem/copy_file.hpp"
#include "utils/perf/timetracer.hpp"

#include "k_range.hpp"
#include "version.hpp"

using fs::make_dir;

namespace spades {
void assemble_genome();
}

struct TimeTracerRAII {
    TimeTracerRAII(llvm::StringRef program_name,
                   unsigned granularity = 500,
                   const std::string &prefix = "", const std::string &suffix = "") {
        time_trace_file_ = prefix + "spades_time_trace_" + suffix + ".json";
        llvm::timeTraceProfilerInitialize(granularity, program_name);
    }
    ~TimeTracerRAII() {
        if (auto E = llvm::timeTraceProfilerWrite(time_trace_file_, "spades-core")) {
            handleAllErrors(std::move(E),
                            [&](const llvm::StringError &SE) {
                                ERROR("" << SE.getMessage() << "\n");
                            });
            return;
        } else {
            INFO("Time trace is written to: " << time_trace_file_);
        }
        llvm::timeTraceProfilerCleanup();
    }

    std::string time_trace_file_;
};

void load_config(const std::vector<std::string>& cfg_fns) {
    for (const auto& s : cfg_fns) {
        fs::CheckFileExistenceFATAL(s);
    }

    cfg::create_instance(cfg_fns);

    make_dir(cfg::get().output_dir);
    make_dir(cfg::get().tmp_dir);

    if (cfg::get().checkpoints != debruijn_graph::config::Checkpoints::None)
        make_dir(cfg::get().output_saves);

    make_dir(cfg::get().temp_bin_reads_path);
}

void create_console_logger(const std::string& dir, std::string log_prop_fn) {
    using namespace logging;

    if (!fs::FileExists(log_prop_fn))
        log_prop_fn = fs::append_path(dir, log_prop_fn);

    logger *lg = create_logger(fs::FileExists(log_prop_fn) ? log_prop_fn : "");
    lg->add_writer(std::make_shared<console_writer>());
    //lg->add_writer(std::make_shared<mutex_writer>(std::make_shared<console_writer>()));
    attach_logger(lg);
}

int main(int argc, char **argv) {
    utils::perf_counter pc;

    const size_t GB = 1 << 30;

    srand(42);
    srandom(42);

    try {
        using namespace debruijn_graph;

        std::string cfg_dir = fs::parent_path(argv[1]);

        std::vector<std::string> cfg_fns;
        for (int i = 1; i < argc; ++i) {
           cfg_fns.push_back(argv[i]);
        }

        // read configuration file (dataset path etc.)
        load_config(cfg_fns);

        create_console_logger(cfg_dir, cfg::get().log_filename);
        for (const auto& cfg_fn : cfg_fns)
            INFO("Loaded config from " << cfg_fn);

        VERIFY(cfg::get().K >= runtime_k::MIN_K && cfg::get().K < runtime_k::MAX_K);
        VERIFY(cfg::get().K % 2 != 0);

        utils::limit_memory(cfg::get().max_memory * GB);

        // assemble it!
        START_BANNER("SPAdes");
        INFO("Maximum k-mer length: " << runtime_k::MAX_K);
        INFO("Assembling dataset (" << cfg::get().dataset_file << ") with K=" << cfg::get().K);
        INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << cfg::get().max_threads);
        std::unique_ptr<TimeTracerRAII> traceraii;
        if (cfg::get().tt.enable || cfg::get().developer_mode) {
            traceraii.reset(new TimeTracerRAII(argv[0],
                                               cfg::get().tt.granularity,
                                               cfg::get().output_dir, std::to_string(cfg::get().K)));
            INFO("Time tracing is enabled");
        }

        TIME_TRACE_SCOPE("spades");
        spades::assemble_genome();
    } catch (std::bad_alloc const &e) {
        std::cerr << "Not enough memory to run SPAdes. " << e.what() << std::endl;
        return EINTR;
    } catch (std::exception const &e) {
        std::cerr << "Exception caught " << e.what() << std::endl;
        return EINTR;
    } catch (...) {
        std::cerr << "Unknown exception caught " << std::endl;
        return EINTR;
    }

    unsigned ms = (unsigned) pc.time_ms();
    unsigned secs = (ms / 1000) % 60;
    unsigned mins = (ms / 1000 / 60) % 60;
    unsigned hours = (ms / 1000 / 60 / 60);
    INFO("Assembling time: " << hours << " hours " << mins << " minutes " << secs << " seconds");

    // OK
    return 0;
}
