//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/core/construction_helper.hpp"

#include "utils/logger/log_writers.hpp"
#include "utils/logger/logger.hpp"
#include "utils/segfault_handler.hpp"

#include "io/graph/gfa_reader.hpp"

#include "version.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <string>

void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

using namespace debruijn_graph;

int main(int argc, char* argv[]) {
    utils::segfault_handler sh;
    utils::perf_counter pc;

    srand(42);
    srandom(42);

    create_console_logger();
    INFO("Starting GFA reader, built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);

    std::string fname(argv[1]);
    ConjugateDeBruijnGraph g(55);

    {
        gfa::GFAReader gfa(fname);
        INFO("Segments: " << gfa.num_edges() << ", links: " << gfa.num_links());

        gfa.to_graph(g);
    }
    
    return 0;
}
