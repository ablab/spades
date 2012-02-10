/*
 * gap_closer.hpp
 *
 *  Created on: Oct 3, 2011
 *
 */

#ifndef GAP_CLOSER_HPP_
#define GAP_CLOSER_HPP_

#include <set>
#include <stack>
#include <unordered_map>

#include "omni/paired_info.hpp"
#include "omni/omni_tools.hpp"
#include "standard.hpp"
#include "omni_labelers.hpp"

namespace debruijn_graph {

const size_t SeedSize = 10;

template<size_t k, class Graph, class SequenceMapper, class Stream>
class GapCloserPairedIndexFiller {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef Seq<SeedSize> SeedMer;

	const Graph &graph_;
	const SequenceMapper& mapper_;
	Stream& stream_;

	unordered_map<EdgeId, pair<EdgeId, int> > OutTipMap;
	unordered_map<EdgeId, pair<EdgeId, int> > InTipMap;
	set<int> InTipsIds;
	set<int> OutTipsIds;
	unordered_map<SeedMer, int> OutHashes;


	size_t CorrectLength(Path<EdgeId> path, size_t idx) {
		size_t answer = graph_.length(path[idx]);
		if (idx == 0)
			answer -= path.start_pos();
		if (idx == path.size() - 1)
			answer -= graph_.length(path[idx]) - path.end_pos();
		return answer;
	}


	Path<EdgeId> ConvertToPath(MappingPath<EdgeId> mp){
		vector<EdgeId> passed;
		for (size_t i = 0; i < mp.size() ;i++){
			passed.push_back(mp[i].first);
		}
		return Path<EdgeId>(passed, -1, passed.size());
	}
	Path<EdgeId> ConvertToPath(Path<EdgeId> mp){
		return mp;
	}



	void ProcessPairedRead(omnigraph::PairedInfoIndex<Graph> &paired_index,
			const io::PairedRead& p_r) {
		Sequence read1 = p_r.first().sequence();
		Sequence read2 = p_r.second().sequence();

		Path<EdgeId> path1 = ConvertToPath(mapper_.MapSequence(read1));
		Path<EdgeId> path2 = ConvertToPath(mapper_.MapSequence(read2));
//		size_t read_distance = p_r.distance();
		for (size_t i = 0; i < path1.size(); ++i) {
//			pair<EdgeId, MappingRange> mapping_edge_1 = path1[i];
			typename unordered_map<EdgeId, pair<EdgeId, int> >::iterator OutTipIter = OutTipMap.find(path1[i]);
			if (OutTipIter != OutTipMap.end()){
				for (size_t j = 0; j < path2.size(); ++j) {
//					pair<EdgeId, MappingRange> mapping_edge_2 = path2[j];
					typename unordered_map<EdgeId, pair<EdgeId, int> >::iterator InTipIter = InTipMap.find(path2[j]);
					if (InTipIter != InTipMap.end()) {
						paired_index.AddPairInfo(PairInfo<EdgeId>(OutTipIter->second.first,	InTipIter->second.first, 100, 1 , 0.));
					}
				}
			}
		}
	}


	void PrepareSeedHashMaps() {
		for (auto iterator = graph_.SmartEdgeBegin(); !iterator.IsEnd(); ) {
			EdgeId edge = *iterator;

			if (graph_.OutgoungEdgeCount(graph_.EdgeEnd(edge)) == 0){
    			Sequence seq1 = graph_.EdgeNucls(edge);
				SeedMer CurSeed = seq1.Subseq(seq1.size() - SeedSize);
				OutHashes.insert(make_pair(CurSeed, graph_.int_id(edge)));
			}
			++iterator;
		}

	}


	void CheckInTipsByHashMaps() {
			for (auto iterator = graph_.SmartEdgeBegin(); !iterator.IsEnd(); ) {
				EdgeId edge = *iterator;
				if (graph_.IncommingEdgeCount(graph_.EdgeStart(edge)) == 0){
	    			Sequence seq1 = graph_.EdgeNucls(edge);
	    			size_t SeqShift = 0;
	    			while (SeqShift + SeedSize < k){
	    				SeedMer CurSeed = seq1.Subseq(SeqShift, SeqShift + SeedSize);
	    				typename unordered_map<EdgeId, pair<EdgeId, int> >::iterator OutTipIter = OutHashes.find(CurSeed);
	    				if (OutTipIter != OutHashes.end()){

	    				}
	    				SeqShift ++;

	    			}
				}
				++iterator;
			}

		}


	void PrepareShiftMaps() {


		stack<pair<EdgeId, int>> edge_stack;
		for (auto iterator = graph_.SmartEdgeBegin(); !iterator.IsEnd(); ) {
			EdgeId edge = *iterator;

			if (graph_.IncomingEdgeCount(graph_.EdgeStart(edge)) == 0){
//				INFO("In tip check");

				InTipMap.insert(make_pair(edge, make_pair(edge, 0)));
				edge_stack.push(make_pair(edge, 0));

				while (edge_stack.size() > 0){
					pair<EdgeId, int> checking_pair = edge_stack.top();
					edge_stack.pop();

					if (graph_.IncomingEdgeCount(graph_.EdgeEnd(checking_pair.first)) == 1){
						if ( graph_.OutgoingEdges(graph_.EdgeEnd(checking_pair.first)).size()>0){
							vector<EdgeId> vec = graph_.OutgoingEdges(graph_.EdgeEnd(checking_pair.first));
							for(size_t i = 0; i < vec.size(); i++){
								EdgeId Cur_edge = vec[i];
								InTipMap.insert(make_pair(Cur_edge , make_pair(edge, graph_.length(checking_pair.first)+ checking_pair.second)));
								edge_stack.push(make_pair(Cur_edge, graph_.length(checking_pair.first) + checking_pair.second));
	//							INFO("Put in stack <"<<Cur_edge<<", "<<graph_.length(edge)<<">");

							}
						}
					}
				}
//				INFO("In tip check complete");
			}


			if (graph_.OutgoingEdgeCount(graph_.EdgeEnd(edge)) == 0){
//				INFO("Out tip check");
				OutTipMap.insert(make_pair(edge, make_pair(edge, 0)));
				edge_stack.push(make_pair(edge, 0));

				while (edge_stack.size() > 0){
					pair<EdgeId, int> checking_pair = edge_stack.top();
					edge_stack.pop();

					if (graph_.OutgoingEdgeCount(graph_.EdgeStart(checking_pair.first)) == 1){
						if ( graph_.IncomingEdges(graph_.EdgeStart(checking_pair.first)).size()>0){
							vector<EdgeId> vec = graph_.IncomingEdges(graph_.EdgeStart(checking_pair.first));
							for(size_t i = 0; i < vec.size(); i++){
								EdgeId Cur_edge = vec[i];
								OutTipMap.insert(make_pair(Cur_edge, make_pair(edge, graph_.length(Cur_edge) + checking_pair.second)));
								edge_stack.push(make_pair(Cur_edge, graph_.length(Cur_edge) + checking_pair.second));
								//INFO("Put in stack <"<<Cur_edge<<", "<<graph_.length(Cur_edge)<<">");
							}
						}
					}

				}
//				INFO("In tip check complete");

			}
			++iterator;
		}
	}


public:

	GapCloserPairedIndexFiller(const Graph &graph, const SequenceMapper& mapper,
			Stream& stream) :
			graph_(graph), mapper_(mapper), stream_(stream) {

	}

	/**
	 * Method reads paired data from stream, maps it to genome and stores it in this PairInfoIndex.
	 */
	void FillIndex(omnigraph::PairedInfoIndex<Graph> &paired_index) {
		INFO("Before prepare");
		PrepareShiftMaps();
		INFO("after prepare");
		stream_.reset();
		while (!stream_.eof()) {
			io::PairedRead p_r;
			stream_ >> p_r;
			ProcessPairedRead(paired_index, p_r);
		}
	}

};

/*template<class Graph>
void PreparePairInfoForGapCloser(Graph& g, omnigraph::PairedInfoIndex<Graph>& paired_info){
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;








}
//vector<EdgeId> IncomingEdges(VertexId v)
*/
template<class Graph, class SequenceMapper>
void CloseShortGaps(Graph& g, omnigraph::PairedInfoIndex<Graph> paired_info, EdgesPositionHandler<Graph> edges_pos, int MimimalIntersection , const SequenceMapper& mapper){
	INFO("Closing short gaps...");
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef vector<PairInfo<EdgeId>> PairInfos;
	DistanceCounter<Graph> distanceTool(g);
	int gaps_filled = 0;
	int gaps_checked = 0;
    for (auto pi_iter = paired_info.begin(); pi_iter != paired_info.end(); ++pi_iter){
    	PairInfos cur_infos = *pi_iter;
    	if (cur_infos.size() > 0){
            VertexId endOfFirstEdge = g.EdgeEnd(cur_infos[0].first);
            VertexId startOfSecondEdge = g.EdgeStart(cur_infos[0].second);
//            if (!distanceTool.IsReachable(endOfFirstEdge, startOfSecondEdge))
              if (cur_infos[0].first != cur_infos[0].second)
              {
            	bool possible_distance = false;
             	size_t k_ = debruijn_graph::K;
             	int cur_gap;
                for(size_t i = 0; i < cur_infos.size(); ++i){
            	   	cur_gap = cur_infos[i].d - g.length(cur_infos[0].first);
//            	 	if (  cur_gap >= 0 && cur_gap <= int(k_ - MimimalIntersection))
                        if ( cur_infos[i].d == 100) 
          	 		possible_distance = true;

            	}
                if (possible_distance){
    				gaps_checked ++;
        			Sequence seq1 = g.EdgeNucls(cur_infos[0].first);
        			Sequence seq2 = g.EdgeNucls(cur_infos[0].second);
           			int best_lev = 100;
           			for (cur_gap = 0; cur_gap <= int(k_ - MimimalIntersection); cur_gap++){
            			if (seq1.Subseq(seq1.size()- k_ +  cur_gap) == seq2.Subseq(0, k_ -  cur_gap)){
            				best_lev = 0;
//                			INFO("YES");
                  			INFO("possible short gap between "<<g.int_id(cur_infos[0].first)<<" and "<<g.int_id(cur_infos[0].second));
                  			INFO("with positions "<<edges_pos.str(cur_infos[0].first)<<"       "<<edges_pos.str(cur_infos[0].second));
                   			INFO("and sequences "<<seq1.Subseq(seq1.size()- k_).str()<<"  "<<seq2.Subseq(0, k_).str());
            				Sequence edge_sequence = seq1.Subseq(seq1.size()- k_) + seq2.Subseq(k_ -  cur_gap, k_);
                			INFO("Gap filled: Gap size = "<<cur_gap<<"  Result seq "<< edge_sequence.str());
                			MappingPath<EdgeId> path1 = mapper.MapSequence(edge_sequence);
                			if (path1.size() > 0) {
                				WARN("Filled k-mer already present in graph");
                			}
                			else {
                				g.AddEdge(endOfFirstEdge, startOfSecondEdge, edge_sequence);
                				gaps_filled ++;
                			}
        					break;
            			} else {
            				int lev = (int)EditDistance(seq1.Subseq(seq1.size()- k_ +  cur_gap), seq2.Subseq(0, k_ -  cur_gap));
            				if (best_lev > lev) best_lev = lev;
            			}

            		}
           			INFO("Best edit distance is "<<best_lev);
            	}
            }
    	}
    }
    INFO("Closing short gaps complete. Total filled " << gaps_filled<<" gaps after checking "<<gaps_checked);
    omnigraph::Compressor<Graph> compressor(g);
    compressor.CompressAllVertices();
}


template<size_t k>
void FirstCloseGap(PairedReadStream& stream, conj_graph_pack& gp){

	typedef NewExtendedSequenceMapper<k + 1, Graph> SequenceMapper;
	stream.reset();
	INFO("First time closing gaps");
	SequenceMapper mapper(gp.g, gp.index, gp.kmer_mapper);
	GapCloserPairedIndexFiller<k + 1, Graph, SequenceMapper, PairedReadStream> gcpif(gp.g, mapper,
			stream);
	paired_info_index gc_paired_info_index(gp.g);
	INFO("FillIndex");
	gcpif.FillIndex(gc_paired_info_index);
	INFO("Close gaps");
	CloseShortGaps(gp.g, gc_paired_info_index, gp.edge_pos, 10, mapper);

}





}





/*


template<class Graph>
void GlueTips(Graph& g, int MimimalIntersection){
	map<Sequence, EdgeId> Starts;
	map<Sequence, EdgeId> Ends;
	set<EdgeId> BlockedEdges;
	set<Sequence> BlockedKmers;
	for (auto iterator = graph_.SmartEdgeBegin(); !iterator.IsEnd(); ) {
		EdgeId tip = *iterator;
		TRACE("Checking edge for being tip " << tip);
		if (IsStartTip(tip)){
			Sequence firstPart = ExtractFirstKMer(tip, CurrentGlueTipK);
			if (BlockedKMers.find(firstPart) != BlockedKmers.end()){
				BlockedEdges.insert(tip);
			}
			else
			if (Starts.find(firstPart) != Starts.end()){
				BlockedEdges.insert()
			}
		}
		TRACE("Use next edge");
		++iterator;
	}

}
*/
#endif /* GAP_CLOSER_HPP_ */
