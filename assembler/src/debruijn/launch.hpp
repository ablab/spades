//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * launch.hpp
 *
 *  Created on: May 6, 2011
 *      Author: sergey
 */

#pragma once

#include "standard.hpp"
#include "omni/visualization/visualization_utils.hpp"

#include "omni/edge_labels_handler.hpp"
#include "graph_construction.hpp"
#include "graph_simplification.hpp"
#include "repeat_resolving.hpp"
#include "omni/omni_tools.hpp"
#include "omni/id_track_handler.hpp"
#include "omni/edges_position_handler.hpp"

#include "de/paired_info.hpp"
#include "de/distance_estimation.hpp"

#include "debruijn_graph.hpp"
#include "config_struct.hpp"
#include "debruijn_stats.hpp"
#include "graphio.hpp"
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
