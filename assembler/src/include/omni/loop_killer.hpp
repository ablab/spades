//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "splitters.hpp"
namespace omnigraph {

template<class Graph>
class AbstractLoopKiller {
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
private:
	Graph &graph_;

	VertexId FindStart(set<VertexId> component_set) {
		VertexId result;
		for(auto it = component_set.begin(); it != component_set.end(); ++it) {
			vector<EdgeId> incoming = graph_.IncomingEdges(*it);
			for(auto eit = incoming.begin(); eit != incoming.end(); ++eit) {
				if(component_set.count(graph_.EdgeStart(*eit)) == 0) {
					if(result != VertexId()) {
						return VertexId();
					}
					result = *it;
				}
			}
		}
		return result;
	}

	VertexId FindFinish(set<VertexId> component_set) {
		VertexId result;
		for(auto it = component_set.begin(); it != component_set.end(); ++it) {
			for (auto I = graph_.out_begin(*it), E = graph_.out_end(*it); I != E; ++I) {
				if (component_set.count(graph_.EdgeEnd(*I)) == 0) {
					if (result != VertexId()) {
						return VertexId();
					}
					result = *it;
				}
			}
		}
		return result;
	}

protected:
	const size_t splitting_edge_length_;
	const size_t max_component_size_;

	Graph &g() {
		return graph_;
	}

public:

	AbstractLoopKiller(Graph &graph, size_t splitting_edge_length,
			size_t max_component_size) :
			graph_(graph), splitting_edge_length_(splitting_edge_length), max_component_size_(
					max_component_size) {
	}

	virtual ~AbstractLoopKiller() {
	}



	void KillAllLoops() {
	    shared_ptr<GraphSplitter<Graph>> splitter_ptr = LongEdgesExclusiveSplitter<Graph>(graph_, splitting_edge_length_);
	    GraphSplitter<Graph> &splitter = *splitter_ptr;
		while(splitter.HasNext()) {
		    set<VertexId> component_set = splitter.Next().vertices();
			if(component_set.size() > max_component_size_)
				continue;
			VertexId start = FindStart(component_set);
			VertexId finish = FindFinish(component_set);
			if(start == VertexId() || finish == VertexId()) {
				continue;
			}
			KillLoop(start, finish, component_set);
		}
		Compressor<Graph> compressor(graph_);
		compressor.CompressAllVertices();
	}

	virtual void KillLoop(VertexId start, VertexId finish, const set<VertexId> &component) = 0;
};

template<class Graph>
class SimpleLoopKiller : public AbstractLoopKiller<Graph> {

public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

private:
	vector<EdgeId> FindPath(VertexId start, VertexId finish, const set<VertexId> &component) {
		set<VertexId> was;
		vector<EdgeId> rr = FindPath(start, finish, was, component);
		return vector<EdgeId>(rr.rbegin(), rr.rend());
	}

	vector<EdgeId> FindPath(VertexId start, VertexId finish, set<VertexId> &was, const set<VertexId> &component) {
		was.insert(start);
		if (start == finish)
			return {};
    for (auto I = this->g().out_begin(start), E = this->g().out_end(start); I != E; ++I) {
      EdgeId edge = *I;
			VertexId next = this->g().EdgeEnd(edge);
			if (next == finish) {
				return { edge };
			}
			if (was.count(next) == 0 && component.count(next) != 0) {
				vector<EdgeId> result = FindPath(next, finish, was, component);
				if (result.size() > 0) {
					result.push_back(edge);
					return result;
				}
			}
		}
		return {};
	}

	bool CheckNotMuchRemoved(const set<EdgeId> &edges, const set<VertexId> &component) {
		size_t sum = 0;
		for (auto it = component.begin(); it != component.end(); ++it) {
      for (auto I = this->g().out_begin(*it), E = this->g().out_end(*it); I != E; ++I) {
        EdgeId edge = *I;
				if (component.count(this->g().EdgeEnd(edge)) == 1 && edges.count(edge) == 0 ) {
					if (this->g().length(edge) > 500) {
						return false;
					}
					sum += this->g().length(edge);
				}
			}
		}
//		if(sum <= 3000) {
//			cout << sum << endl;
//		}
		return sum <= 3000;
	}

	void RemoveExtraEdges(const set<EdgeId> &edges, const set<VertexId> &component) {
		vector<VertexId> comp(component.begin(), component.end());
		vector<EdgeId> to_delete;
		for (auto it = comp.begin(); it != comp.end(); ++it) {
      for (auto I = this->g().out_begin(*it), E = this->g().out_end(*it); I != E; ++I) {
        EdgeId edge = *I;
				if (component.count(this->g().EdgeEnd(edge)) == 1 && edges.count(edge) == 0) {
					to_delete.push_back(edge);
				}
			}
		}
    
		SmartSetIterator<Graph, EdgeId> s(this->g(), to_delete.begin(), to_delete.end());
		while (!s.IsEnd()) {
			this->g().DeleteEdge(*s);
			++s;
		}
	}

	void RemoveIsolatedVertices(set<VertexId> component) {
		SmartSetIterator<Graph, VertexId> s(this->g(), component.begin(), component.end());
		while (!s.IsEnd()) {
			if (this->g().IsDeadStart(*s) && this->g().IsDeadEnd(*s)) {
				this->g().DeleteVertex(*s);
			}
			++s;
		}
	}

	bool CheckStrong(const set<VertexId> &component) {
		VertexId v = *(component.begin());
		for(auto it = component.begin(); it != component.end(); ++it) {
			if(v != *it && (FindPath(v, *it, component).size() == 0 || FindPath(*it, v, component).size() == 0)) {
				return false;
			}
		}
		return true;
	}

public:
	SimpleLoopKiller(Graph &graph, size_t splitting_edge_length, size_t max_component_size) :
			AbstractLoopKiller<Graph>(graph, splitting_edge_length, max_component_size) {
	}

	virtual void KillLoop(VertexId start, VertexId finish, const set<VertexId> &component) {
		vector<EdgeId> path = FindPath(start, finish, component);
		set<EdgeId> edges(path.begin(), path.end());
		if(path.size() > 0 || start == finish) {
//			if(start != finish || component.size() > 2)
			if(/*!CheckStrong(component) || */!CheckNotMuchRemoved(edges, component)) {
				return;
			}
/*
			cout << this->g().int_id(start) << " " << this->g().int_id(finish) << endl;
			cout << this->g().VertexNucls(start) << endl;

			for(auto it = component.begin(); it != component.end(); ++it) {
				vector<EdgeId> outgoing = this->g().OutgoingEdges(*it);
				for(auto eit = outgoing.begin(); eit != outgoing.end(); ++eit) {
					if(component.count(this->g().EdgeEnd(*eit)) == 1) {
						cout << this->g().int_id(*it) << " -> " << this->g().int_id(this->g().EdgeEnd(*eit)) << " : " << this->g().length(*eit) << " : " << edges.count(*eit) << " : " << this->g().int_id(*eit) << endl;
					}
				}
			}
*/
			RemoveExtraEdges(edges, component);
			RemoveIsolatedVertices(component);
		}
	}

};
}
