/*
 * long_read_mapper.hpp
 *
 *  Created on: Jun 17, 2013
 *      Author: andrey
 */

#ifndef LONG_READ_MAPPER_HPP_
#define LONG_READ_MAPPER_HPP_

#include "pe_utils.hpp"

namespace path_extend {

class SimpleLongReadMapper {

    conj_graph_pack& gp_;

    ExtendedSequenceMapper<Graph> mapper_;

    SameEdgeDeletionCorrector same_edge_corr_;

    DijkstraSearcher ps_;

    CloseGapsCorrector gap_closer_;

public:
    SimpleLongReadMapper(conj_graph_pack& conj_gp): gp_(conj_gp),
        mapper_(gp_.g, gp_.index, gp_.kmer_mapper, gp_.k_value + 1),
        same_edge_corr_(gp_.g),
        ps_(gp_.g),
        gap_closer_(gp_.g, &ps_)
    {
        paths_searcher_config conf;
        conf.depth_neigh_search = 5; // max path len (in edges)
        conf.max_len_path = 100000;  // max path len (in k-mers)
        conf.max_num_vertices = 100; // max number of visited vertices
        ps_.Initialize(conf);
    }

    template<class SingleRead>
    void ProcessLib(io::IReader<SingleRead>& stream, PathStorage<Graph>& storage) {
    	while (!stream.eof()) {
            SingleRead r;
            stream >> r;

            MappingPath<EdgeId> path;
            path.join(mapper_.MapSequence(r.sequence()));

            SimpleMappingContig mc(r.sequence(), path);
            MappingContig * dc = same_edge_corr_.Correct(&mc);

            MappingContig * gc = gap_closer_.Correct(dc);

            storage.AddPath(gc->PathSeq(), 1, true);
        }
    }

};


}


#endif /* LONG_READ_MAPPER_HPP_ */
