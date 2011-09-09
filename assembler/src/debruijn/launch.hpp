/*
 * launch.hpp
 *
 *  Created on: May 6, 2011
 *      Author: sergey
 */

#pragma once

#include "standart.hpp"
#include "omni/visualization_utils.hpp"

//#include "debruijn_graph.hpp"
#include "omni/edge_labels_handler.hpp"
#include "omni/paired_info.hpp"
#include "graph_construction.hpp"
#include "graph_simplification.hpp"
#include "coverage_handler.hpp"
#include "repeat_resolving.hpp"
#include "omni/omni_tools.hpp"
#include "seq_map.hpp"
#include "omni/ID_track_handler.hpp"
#include "omni/edges_position_handler.hpp"

#include "new_debruijn.hpp"
#include "config_struct.hpp"
#include "debruijn_stats.hpp"
#include "graphio.hpp"
#include "rectangleRepeatResolver.hpp"
#include "omni/distance_estimation.hpp"
#include "omni/loop_resolver.hpp"
#include "check_tools.hpp"
#include "copy_file.hpp"

#include "graph_pack.hpp"
#include "repeat_resolving_routine.hpp"

#include "omni_labelers.hpp"

namespace debruijn_graph
{

void assembly_genome(PairedReadStream& stream, const Sequence& genome)
{
    INFO("Assembly Genome Started");
    TRACE("Starting from: " << debruijn_config::working_stage_name(cfg::get().entry_point));

    make_repeat_resolving(stream, genome);

    INFO("Assembly Genome Finished");
}

}
