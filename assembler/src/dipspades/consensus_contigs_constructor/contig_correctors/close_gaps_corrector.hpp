#pragma once

#include "abstract_contig_corrector.hpp"

using namespace debruijn_graph;

namespace dipspades {

class CloseGapsCorrector : public AbstractContigCorrector{

	vector<int> incorr_contigs;
	size_t num_corr;

	vector<EdgeId> ClosePathGap(vector<EdgeId> path, vector<size_t> gap_index){
		vector<EdgeId> new_path;
		size_t current_gap = 0;
		for(size_t i = 0; i < path.size() - 1; i++){
			EdgeId cur_edge = path[i];
			new_path.push_back(cur_edge);
			if(i == gap_index[current_gap]){
				VertexId start = g_.EdgeEnd(cur_edge);
				VertexId end = g_.EdgeStart(path[i + 1]);
				auto dijkstra = DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, dsp_cfg::get().pbr.max_bulge_nucls_len);
				dijkstra.run(start);
				if(dijkstra.DistanceCounted(end)){
					vector<EdgeId> add_path = dijkstra.GetShortestPathTo(end);
					for(auto e = add_path.begin(); e != add_path.end(); e++)
						if(g_.EdgeStart(*e) != g_.EdgeEnd(*e))
							new_path.push_back(*e);
				}
				else{
					// second attempt
					VertexId prev_start = g_.EdgeStart(cur_edge);
					dijkstra.run(prev_start);
	               if(dijkstra.DistanceCounted(end)){
	                    vector<EdgeId> add_path = dijkstra.GetShortestPathTo(end);
						new_path.erase(new_path.begin() + new_path.size() - 1);
                	    for(auto e = add_path.begin(); e != add_path.end(); e++)
                	    	if(g_.EdgeStart(*e) != g_.EdgeEnd(*e))
                	    		new_path.push_back(*e);
					}
				}
				current_gap++;
			}
		}
		new_path.push_back(path[path.size() - 1]);
		return new_path;
	}

	size_t CountContigsWithGaps(ContigStoragePtr storage) {
		size_t contigs_with_gaps = 0;
		for(size_t i = 0; i < storage->Size(); i++)
			if(!IsPathConnected(g_, (*storage)[i]->path_seq()))
				contigs_with_gaps++;
		return contigs_with_gaps;
	}

	void ProcessContigs(ContigStoragePtr storage) {
		double processed_perc = 0.1;
		double step = 0.1;
		for(size_t i = 0; i < storage->Size(); i++) {
			storage->ReplaceContig(Correct((*storage)[i]), i);
			double cur_process_perc = static_cast<double>(i) / static_cast<double>(storage->Size());
			if(cur_process_perc > processed_perc) {
				while(processed_perc + step <= cur_process_perc)
					processed_perc += step;
				INFO(ToString(processed_perc * 100.0) << "% contigs were processed");
				processed_perc += step;
			}
		}
		INFO("100% contigs were processed");
	}

public:
	CloseGapsCorrector(Graph &g) : AbstractContigCorrector(g) {
		num_corr = 0;
	}

	virtual ContigStoragePtr Correct(ContigStoragePtr storage){

		INFO(ToString(CountContigsWithGaps(storage)) << " contigs from " <<
				ToString(storage->Size()) << " have gaps before correction");

		ProcessContigs(storage);

		INFO(ToString(num_corr) + " contigs from " + ToString(storage->Size()) + " were corrected");
		INFO(ToString(storage->Size() - num_corr) + " contigs from " << storage->Size() << " have gaps after correction");
		return storage;
	}

	virtual MappingContigPtr Correct(MappingContigPtr contig){
		vector<EdgeId> path = contig->path_seq();
		if(path.size() <= 1){
			num_corr++;
			return contig;
		}
		vector<size_t> gap_indexes;
		for(size_t i = 0; i < path.size() - 1; i++){
			EdgeId e1 = path[i];
			EdgeId e2 = path[i + 1];
			if(!AreEdgesConnected(g_, e1, e2)){
				gap_indexes.push_back(i);
			}
		}

		if(gap_indexes.size() != 0){
			vector<EdgeId> new_path = ClosePathGap(path, gap_indexes);
			if(IsPathConnected(g_, new_path))
				num_corr++;
			return MappingContigPtr(new ReplacedPathMappingContig(contig, new_path));
		}
		num_corr++;
		TRACE("Contig " << contig->id() << " has a gap");
		TRACE("Contig path: " << SimplePathWithVerticesToString(g_, contig->path_seq()));
		return contig;
	}
};

}
