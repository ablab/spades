/*
 * launch.hpp
 *
 *  Created on: May 6, 2011
 *      Author: sergey
 */

#pragma once

#include "standard.hpp"
#include "omni/visualization_utils.hpp"

//#include "debruijn_graph.hpp"
#include "omni/edge_labels_handler.hpp"
#include "omni/paired_info.hpp"
#include "graph_construction.hpp"
#include "graph_simplification.hpp"
#include "repeat_resolving.hpp"
#include "omni/omni_tools.hpp"
#include "seq_map.hpp"
#include "omni/id_track_handler.hpp"
#include "omni/edges_position_handler.hpp"

#include "new_debruijn.hpp"
#include "config_struct.hpp"
#include "debruijn_stats.hpp"
#include "graphio.hpp"
#include "omni/distance_estimation.hpp"
//#include "omni/advanced_distance_estimation.hpp"
#include "omni/loop_resolver.hpp"
#include "check_tools.hpp"
#include "copy_file.hpp"

#include "graph_pack.hpp"
#include "repeat_resolving_routine.hpp"

#include "omni_labelers.hpp"

namespace debruijn_graph
{

void assemble_genome()
{
    INFO("SPAdes started");
    INFO("Starting from stage: " << debruijn_config::working_stage_name(cfg::get().entry_point));

    exec_repeat_resolving();

    INFO("SPAdes finished");

}
}
