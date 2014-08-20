#pragma once

#include "../utils/bulge_utils.hpp"

using namespace debruijn_graph;

namespace dipspades {

typedef shared_ptr<BaseBulge> BulgePtr;

class BulgeChain {
	const Graph &graph_;

	vector<EdgeId> conservative_regions_;
	vector<BulgePtr> bulges_;

	enum chain_elem_type { bulge_region, edge_region };

	struct chain_element 	{
		size_t index;
		chain_elem_type elem_type;

		chain_element(chain_elem_type new_elem) :
			index(-1),
			elem_type(new_elem) { }
	};

	vector<chain_element> chain_elems_;

public:
	BulgeChain(const Graph &graph) :
		graph_(graph) { }

	size_t EdgesNumber() { return conservative_regions_.size() + bulges_.size() * 2; }

	double PercentIdentity() {
		double sum_diversity = 0;
		size_t sum_length = 0;
		for (auto b = bulges_.begin(); b != bulges_.end(); b++) {
			sum_diversity += (*b)->relative_align() * double((*b)->BulgeLength());
			sum_length += (*b)->BulgeLength();
		}

		for (auto e = conservative_regions_.begin(); e != conservative_regions_.end(); e++) {
			sum_diversity += double(graph_.length(*e));
			sum_length += graph_.length(*e);
		}

		return sum_diversity / double(sum_length);
	}

	size_t Length() {
		size_t length = 0;
		for (auto b = bulges_.begin(); b != bulges_.end(); b++) {
			length += (*b)->BulgeLength();
		}

		for (auto e = conservative_regions_.begin(); e != conservative_regions_.end(); e++) {
			length += graph_.length(*e);
		}

		return length;
	}

	void AddBulgeToStart(BulgePtr bulge) {
		bulges_.insert(bulges_.begin(), bulge);
		chain_elems_.insert(chain_elems_.begin(), chain_element(bulge_region));
	}

	void AddBulgeToEnd(BulgePtr bulge) {
		bulges_.push_back(bulge);
		chain_elems_.push_back(chain_element(bulge_region));
	}

	void AddEdgeToStart(EdgeId edge) {
		conservative_regions_.insert(conservative_regions_.begin(), edge);
		chain_elems_.insert(chain_elems_.begin(), chain_element(edge_region));
	}

	void AddEdgeToEnd(EdgeId edge) {
		conservative_regions_.push_back(edge);
		chain_elems_.push_back(chain_element(edge_region));
	}

	typedef vector<BulgePtr>::iterator bulge_iterator;

	bulge_iterator BulgeBegin() { return bulges_.begin(); }

	bulge_iterator BulgeEnd() { return bulges_.end(); }

	string ToString() {
		size_t bulge_ptr = 0;
		size_t edge_ptr = 0;
		stringstream sstream;

		sstream << "Chain length: " << Length() << endl;
		sstream << "% identity: " << PercentIdentity() << endl;

		for (auto elem = chain_elems_.begin(); elem != chain_elems_.end(); elem++) {
			if (elem->elem_type == bulge_region) {
				sstream << "Bulge: " << bulges_[bulge_ptr]->BulgeToString() << endl;
				bulge_ptr++;
			}
			else {
				sstream << "Conservative region: " << graph_.int_id(conservative_regions_[edge_ptr]) << endl;
				edge_ptr++;
			}
		}
		return sstream.str();
	}
};

class BugleChainFinder {
	size_t min_length_;
	const Graph &graph_;

	double rel_length_threshold_;
	double rel_seq_threshold_;

	set<string> processed_bugles_;

	bool BulgeIsGood(BulgePtr bulge) {
		string str_id = bulge->StrId();
		return processed_bugles_.find(str_id) == processed_bugles_.end() and bulge->IsBulgeDiploid(rel_length_threshold_, rel_seq_threshold_); // and BulgeIsDiploid()
	}

	void ForwardBulgeProcessor(BulgePtr seed_bulge, BulgeChain &bulge_chain) {
		bool chain_is_active = true;
		VertexId cur_vertex = seed_bulge->end_vertex();
		while(chain_is_active) {
			auto outgoing_edges = graph_.OutgoingEdges(cur_vertex);

			// conservative region is outgoing
			if(graph_.OutgoingEdgeCount(cur_vertex) == 1) {
				EdgeId cons_region = *(outgoing_edges.begin());
				bulge_chain.AddEdgeToEnd(cons_region);
				cur_vertex = graph_.EdgeEnd(cons_region);
			}
			// > 1 outgoing edges
			else {
				if(graph_.OutgoingEdgeCount(cur_vertex) == 2) {
					EdgeId edge1 = *(outgoing_edges.begin());
					EdgeId edge2 = *(++outgoing_edges.begin());

					if(graph_.EdgeEnd(edge1) != graph_.EdgeEnd(edge2))
						break; // edges do not form a bulge, chain is finished

					auto bulge = shared_ptr<BaseBulge>(new Bulge(graph_, graph_.k(), edge1, edge2));
					if (bulge->IsBulgeDiploid(rel_length_threshold_, rel_seq_threshold_))
						bulge_chain.AddBulgeToEnd(bulge);
				}
				else
					break; // chain is finished
			}
		}
	}

	void BackwardBulgeProcessor(BulgePtr seed_bulge, BulgeChain &bulge_chain) {
		bool chain_is_active = true;
		VertexId cur_vertex = seed_bulge->start_vertex();
		while(chain_is_active) {
			auto incoming_edges = graph_.IncomingEdges(cur_vertex);

			// conservative region is incoming
			if(graph_.IncomingEdgeCount(cur_vertex) == 1) {
				EdgeId cons_region = *(incoming_edges.begin());
				bulge_chain.AddEdgeToStart(cons_region);
				cur_vertex = graph_.EdgeStart(cons_region);
			}

			// > 1 incoming edges
			else {
				if(graph_.IncomingEdgeCount(cur_vertex) == 2) {
					EdgeId edge1 = *(incoming_edges.begin());
					EdgeId edge2 = *(++incoming_edges.begin());

					if(graph_.EdgeStart(edge1) != graph_.EdgeStart(edge2))
						break; // edges do not form a bulge, chain is finished

					auto bulge = shared_ptr<BaseBulge>(new Bulge(graph_, graph_.k(), edge1, edge2));
					if (bulge->IsBulgeDiploid(rel_length_threshold_, rel_seq_threshold_))
						bulge_chain.AddBulgeToStart(bulge);
				}
				else
					break; // chain is finished
				}
		}
	}

	BulgeChain FindChainByBulge(BulgePtr seed_bulge) {
		BulgeChain bulge_chain(graph_);
		ForwardBulgeProcessor(seed_bulge, bulge_chain);
		BackwardBulgeProcessor(seed_bulge, bulge_chain);
		return bulge_chain;
	}

public:
	BugleChainFinder(const Graph &graph,
			size_t min_length,
			double rel_len_threshold,
			double rel_seq_threshold) :
		min_length_(min_length),
		graph_(graph),
		rel_length_threshold_(rel_len_threshold),
		rel_seq_threshold_(rel_seq_threshold){ }

	void FindGoodChains() {
		for(auto edge = SmartEdgeIterator<Graph>(graph_); !edge.IsEnd(); ++edge) {
			vector<EdgeId> edges = graph_.GetEdgesBetween(graph_.EdgeStart(*edge),
					graph_.EdgeEnd(*edge));
				if(edges.size() >= 2){
					auto bulge = shared_ptr<BaseBulge>(new Bulge(graph_, graph_.k(), edges[0], edges[1]));
					string str_id = bulge->StrId();
					TRACE("BulgeId: " << str_id);
					if(BulgeIsGood(bulge)) {
						BulgeChain chain = FindChainByBulge(bulge);
						TRACE("Chain was found");
						for (auto b = chain.BulgeBegin(); b != chain.BulgeEnd(); b++)
							processed_bugles_.insert((*b)->StrId());
						TRACE(chain.ToString());
					}
				}
		}
	}

private:
	DECL_LOGGER("BulgeChainFinder");
};

class HeterozygosityEstimator {
	const conj_graph_pack &graph_pack_;
	BugleChainFinder chain_finder_;

public:
	HeterozygosityEstimator(conj_graph_pack &graph_pack, double rel_len_threshold, double rel_seq_threshold) :
		graph_pack_(graph_pack),
		chain_finder_(graph_pack.g, 5000, rel_len_threshold, rel_seq_threshold) { }

	void Estimate() {
		chain_finder_.FindGoodChains();
	}
};

}
