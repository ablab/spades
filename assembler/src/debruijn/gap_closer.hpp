/*
 * gap_closer.hpp
 *
 *  Created on: Oct 3, 2011
 *
 */

#ifndef GAP_CLOSER_HPP_
#define GAP_CLOSER_HPP_

#include "omni/paired_info.hpp"


template<class Graph>
void CloseShortGaps(Graph& g, omnigraph::PairedInfoIndex<Graph> paired_info, EdgesPositionHandler<Graph> edges_pos){
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef vector<PairInfo<EdgeId>> PairInfos;

	DistanceCounter<Graph> distanceTool(g);
	int gaps_filled = 0;
    for (auto pi_iter = paired_info.begin(); pi_iter != paired_info.end(); ++pi_iter){
    	PairInfos cur_infos = *pi_iter;
    	if (cur_infos.size() > 0){
            VertexId endOfFirstEdge = g.EdgeEnd(cur_infos[0].first);
            VertexId startOfSecondEdge = g.EdgeStart(cur_infos[0].second);
            if (!distanceTool.IsReachable(endOfFirstEdge, startOfSecondEdge))
            {
            	for(size_t i = 0; i < cur_infos.size(); ++i){
                	int cur_gap = cur_infos[i].d - g.length(cur_infos[0].first);
            		if (  cur_gap >= 0 && cur_gap <= 45){
            			Sequence seq1 = g.EdgeNucls(cur_infos[0].first);
            			Sequence seq2 = g.EdgeNucls(cur_infos[0].second);
            			size_t k_ = debruijn_graph::K;
//            			INFO("possible short gap "<<edges_pos.str(cur_infos[0].first)<<"       "<<edges_pos.str(cur_infos[0].second));
//            			INFO("possible short gap "<<seq1.Subseq(seq1.size()- k_ + cur_gap).str()<<"  "<<seq2.Subseq(0, k_ -  cur_gap).str());
            			if (seq1.Subseq(seq1.size()- k_ +  cur_gap) == seq2.Subseq(0, k_ -  cur_gap)){
//                			INFO("YES");
            				Sequence edge_sequence = seq1.Subseq(seq1.size()- k_) + seq2.Subseq(k_ -  cur_gap, k_);
            				INFO("possible short gap "<<edges_pos.str(cur_infos[0].first)<<"       "<<edges_pos.str(cur_infos[0].second));
//                			INFO("possible short gap "<<seq1.Subseq(seq1.size()- k_ + cur_gap).str()<<"  "<<seq2.Subseq(0, k_ -  cur_gap).str());
                			INFO("Gap filled: Gap size = "<<cur_gap<<"  Result seq "<< edge_sequence.str());
            				g.AddEdge(endOfFirstEdge, startOfSecondEdge, edge_sequence);
            				gaps_filled ++;
        					break;
            			}

            		}
            	}
            }
    	}
    }
    INFO("Total filled " << gaps_filled<<" gaps");
    omnigraph::Compressor<Graph> compressor(g);
    compressor.CompressAllVertices();

}



#endif /* GAP_CLOSER_HPP_ */
