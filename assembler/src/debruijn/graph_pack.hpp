#pragma once

#include "omni/edge_labels_handler.hpp"
#include "omni/paired_info.hpp"
#include "omni/ID_track_handler.hpp"
#include "omni/edges_position_handler.hpp"
#include "new_debruijn.hpp"
#include "config_struct.hpp"
#include "graphio.hpp"

namespace debruijn_graph
{

typedef PairedInfoIndex<Graph>  paired_info_index;

struct conj_graph_pack
{
    typedef ConjugateDeBruijnGraph graph_t;

    graph_t                                     g;
	EdgeIndex<debruijn_graph::K + 1, graph_t>   index;
	IdTrackHandler<graph_t>                     int_ids;
	EdgesPositionHandler<graph_t>               edge_pos;
	PairedInfoIndex<graph_t>                    etalon_paired_index;
	KmerMapper<debruijn_graph::K + 1, graph_t>  kmer_mapper;

	Sequence const& genome;

	conj_graph_pack(Sequence const& genome)
	    : g         (debruijn_graph::K)
	    , index     (g)
	    , int_ids   (g)
	    , edge_pos  (g)
	    , etalon_paired_index(g, 0)
	    , kmer_mapper(g)
	    , genome    (genome)
	{
	}
};

struct nonconj_graph_pack
{
    typedef NonconjugateDeBruijnGraph graph_t;

    graph_t                                     g;
    IdTrackHandler<graph_t>                     int_ids;
    EdgesPositionHandler<graph_t>               edge_pos;
    PairedInfoIndex<NCGraph>                    clustered_index;

    nonconj_graph_pack()
        : g                 (K)
        , int_ids           (g)
        , edge_pos          (g)
        , clustered_index   (g)
    {
    }

    nonconj_graph_pack(conj_graph_pack const& gp)
        : g         (K)
        , int_ids   (g)
        , edge_pos  (g)
    {
        fs::path p = cfg::get().output_root / "temp_conversion" / "conj_graph";

        printGraph(gp.g, gp.int_ids, p.string(), paired_info_index(),
                   gp.edge_pos, &gp.etalon_paired_index, &clustered_index);

        scanNCGraph<graph_t>(g, int_ids, p.string(), 0, edge_pos, 0, &clustered_index);
    }
};

} // namespace debruijn_graph
/*
 * graph_pack.hpp
 *
 *  Created on: Aug 18, 2011
 *      Author: sergey
 */

#pragma once

#include "omni/edge_labels_handler.hpp"
#include "omni/paired_info.hpp"
#include "omni/ID_track_handler.hpp"
#include "omni/edges_position_handler.hpp"
#include "new_debruijn.hpp"
#include "config_struct.hpp"
#include "graphio.hpp"

namespace debruijn_graph
{

typedef PairedInfoIndex<Graph>  paired_info_index;

struct conj_graph_pack
{
    typedef ConjugateDeBruijnGraph graph_t;

    graph_t                                     g;
	EdgeIndex<debruijn_graph::K + 1, graph_t>   index;
	IdTrackHandler<graph_t>                     int_ids;
	EdgesPositionHandler<graph_t>               edge_pos;
	PairedInfoIndex<graph_t>                    etalon_paired_index;
	KmerMapper<debruijn_graph::K + 1, graph_t>  kmer_mapper;

	Sequence const& genome;

	conj_graph_pack(Sequence const& genome)
	    : g         (debruijn_graph::K)
	    , index     (g)
	    , int_ids   (g)
	    , edge_pos  (g)
	    , etalon_paired_index(g, 0)
	    , kmer_mapper(g)
	    , genome    (genome)
	{
	}
};

struct nonconj_graph_pack
{
    typedef NonconjugateDeBruijnGraph graph_t;

    graph_t                                     g;
    IdTrackHandler<graph_t>                     int_ids;
    EdgesPositionHandler<graph_t>               edge_pos;
    PairedInfoIndex<NCGraph>                    clustered_index;

    nonconj_graph_pack()
        : g                 (K)
        , int_ids           (g)
        , edge_pos          (g)
        , clustered_index   (g)
    {
    }

    nonconj_graph_pack(conj_graph_pack const& gp)
        : g         (K)
        , int_ids   (g)
        , edge_pos  (g)
    {
        fs::path p = cfg::get().output_root / "temp_conversion" / "conj_graph";
        printGraph(gp.g, gp.int_ids, p.string(), paired_info_index(),
                   gp.edge_pos, &gp.etalon_paired_index, &clustered_index);

        scanNCGraph<graph_t>(g, int_ids, p.string(), 0, edge_pos, 0, &clustered_index);


    }
   /* void nonconj_graph_copy (conj_graph_pack const& cg, nonconj_graph_pack& ncg) {
    	for(auto iter = cg.g.SmartVertexBegin(); !iter.IsEnd(); ++iter) {
    		ncg.g.AddVertex();
    		ncg.
    	}


    }
*/
};

} // namespace debruijn_graph


__________________________
