//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

// CGCE stands for 'colored graph correction evaluation'

#include "compare_standard.hpp"
#include "coloring.hpp"
#include "utils/logger/logger.hpp"
#include "utils/logger/log_writers.hpp"

#include "coordinates_handler.hpp" // These lines are needed for
#include "pipeline/config_struct.hpp"       // colored_graph_construction.hpp
using namespace debruijn_graph;    // Textual inclusion is awesome!

#include "colored_graph_construction.hpp"

using namespace cap;
using namespace logging;

int main(int argc, char** argv) {
    if (argc < 3)
        return -1;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);

    const size_t k = 55;
    std::string uncorrected_fn(argv[1]);
    std::string corrected_fn(argv[2]);

    typedef conj_graph_pack::graph_t Graph;
    typedef conj_graph_pack::index_t Index;

    conj_graph_pack gp(k, ".", 0);

    vector<ContigStream*> streams;
    streams.push_back(new io::FileReadStream(uncorrected_fn));
    streams.push_back(new io::FileReadStream(corrected_fn));
        
    vector<ContigStream*> rc_contigs;
    for (auto it = streams.begin(); it != streams.end(); ++it) {
      rc_contigs.push_back(new RCWrapper(**it));
      rc_contigs.back()->reset();
    }

    ContigStreams rc_read_stream_vector(rc_contigs, true);

    ColorHandler<Graph> coloring(gp.g, streams.size());

    debruijn_config::construction params;
    params.con_mode = construction_mode::con_extention;
    params.early_tc.enable = false;
    params.keep_perfect_loops = true;
    ConstructGraph(k, params, omp_get_max_threads(), rc_read_stream_vector, gp.g, gp.index);
    SplitAndColorGraph(gp, coloring, rc_read_stream_vector);

    std::ofstream red("red.fa");
    std::ofstream blue("blue.fa");
    
    size_t n_red = 0, n_blue = 0;
    
    for (auto I = gp.g.begin(), E = gp.g.end(); I != E; ++I) {
        for (auto eit = gp.g.out_begin(*I); eit != gp.g.out_end(*I); ++eit) {
            if (coloring.Color(*eit) == kRedColorSet)
                red << ">RED_" << ++n_red << "\n"
                    << gp.g.EdgeNucls(*eit).str() << std::endl;
            else if (coloring.Color(*eit) == kBlueColorSet)
                blue << ">BLUE_" << ++n_blue << "\n"
                    << gp.g.EdgeNucls(*eit).str() << std::endl;
        }
    }
}
