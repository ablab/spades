/**
 * condensed_graph.cpp
 *
 *  Created on: Feb 21, 2011
 *      Author: sergey
 */
#include "condensed_graph.hpp"
#include "logging.hpp"
#include "graphVisualizer.hpp"
#include <set>
#include <iostream>
#include <tr1/unordered_map>

using namespace std;

namespace condensed_graph {

void CondensedGraph::FixIncomingOnSplit(Vertex* v, Vertex* v1, Vertex* v2) {
	vector<Vertex*> anc = LeftNeighbours(v);
	for (size_t i = 0; i < anc.size(); ++i) {
		Vertex* ancestor = anc[i];
		if (ancestor == v->complement()) {
			LinkVertices(v1 -> complement(), v1);
		} else if (ancestor == v) {
			LinkVertices(v2, v1);
		} else {
			//trivial case
			LinkVertices(ancestor, v1);
		}
	}
}

void CondensedGraph::FixIncomingOnMerge(Vertex* v1, Vertex* v2, Vertex* v) {
	vector<Vertex*> anc = LeftNeighbours(v1);
	for (size_t i = 0; i < anc.size(); ++i) {
		Vertex* ancestor = anc[i];
		if (ancestor == v1->complement()) {
			LinkVertices(v->complement(), v);
		} else if (ancestor == v2) {
			LinkVertices(v, v);
		} else {
			//trivial case
			LinkVertices(ancestor, v);
		}
	}
}

bool CondensedGraph::CheckIfNoIncoming(Vertex* v) const {
	//this code executes while graph contract is temporarily broken
	//vertex from Anc(v) doesn't necessarily have link to v
	vector<Vertex*> anc = LeftNeighbours(v);
	for (size_t i = 0; i < anc.size(); ++i) {
		Vertex* ancestor = anc[i];
		if (ancestor != v && ancestor != v->complement()) {
			for (size_t j = 0; j < 4; ++j) {
				if (ancestor->right_neighbour(j) == v) {
					return false;
				}
			}
		}
	}
	return true;
}

bool CondensedGraph::CanBeDeleted(Vertex* v) const {
	return CheckIfNoIncoming(v) && CheckIfNoIncoming(v->complement());
}

vector<Vertex*> CondensedGraph::LeftNeighbours(const Vertex* v) const {
	vector<Vertex*> ans;
	Vertex* complement = v->complement();
	for (char i = 3; i >= 0; --i) {
		if (complement->right_neighbour(i) != NULL) {
			ans.push_back(complement->right_neighbour(i)->complement());
		}
	}
	return ans;
}

vector<Vertex*> CondensedGraph::RightNeighbours(const Vertex* v) const {
	vector<Vertex*> ans;
	for (char i = 0; i < 4; ++i) {
		Vertex* right_neighbour = v->right_neighbour(i);
		if (right_neighbour != NULL) {
			ans.push_back(right_neighbour);
		}
	}
	return ans;
}

Vertex* CondensedGraph::AddVertex(const Sequence &nucls) {
	DEBUG("Adding vertex for sequence '" << nucls.str() << "' and its complement '" << (!nucls).str() << "'");

	Vertex* v1 = new Vertex(nucls);
	Vertex* v2 = new Vertex(!nucls);
	v1->set_complement(v2);
	v2->set_complement(v1);
	vertices_.insert(v1);
	vertices_.insert(v2);
	for (size_t i = 0; i < action_handler_list_.size(); ++i) {
		action_handler_list_[i]->HandleAdd(v1);
	}
	return v1;
}

void CondensedGraph::DeleteVertex(Vertex* v) {
	DEBUG("Deleting vertex '" << v->nucls().str() << "' and its complement '" << v->complement()->nucls().str() << "'")

	assert(CanBeDeleted(v));

	Vertex* complement = v->complement();
	vertices_.erase(v);
	vertices_.erase(complement);

	for (size_t i = 0; i < action_handler_list_.size(); ++i) {
		action_handler_list_[i]->HandleDelete(v);
	}

	delete v;
	delete complement;
}

Vertex* CondensedGraph::SplitVertex(Vertex* v, size_t pos) {
	DEBUG("Splitting vertex '" << v->nucls().str() <<"' of size " << v->size() << " at position "<< pos);
	assert(pos <= v->size());

	if (pos == v->size()) {
		return v;
	};

	Sequence nucls = v->nucls();

	Vertex* v1 = AddVertex(nucls.Subseq(0, pos));
	Vertex* v2 = AddVertex(nucls.Subseq(pos - (k_ - 1), nucls.size())); // nucls.size() can be omitted here

	LinkVertices(v1, v2);

	FixIncomingOnSplit(v, v1, v2);

	FixIncomingOnSplit(v->complement(), v2->complement(), v1->complement());

//	action_handler_->HandleSplit(v, pos, v1, v2);

	DeleteVertex(v);

	return v1;
}

Vertex* CondensedGraph::Merge(Vertex* v1, Vertex* v2) {
	DEBUG("Merging vertices '" << v1->nucls().str() << "' and '" << v2->nucls().str() << "' and their complement")
	assert(IsMergePossible(v1, v2));

	Vertex* v = AddVertex(v1->nucls() + v2->nucls().Subseq(k_ - 1));
	FixIncomingOnMerge(v1, v2, v);
	FixIncomingOnMerge(v2->complement(), v1->complement(), v->complement());

//	action_handler_->HandleMerge(v1, v2, v);

	DeleteVertex(v1);
	DeleteVertex(v2);
	return v;
}

void CondensedGraph::AddRightNeighbour(Vertex* v1, Vertex* v2) {
	v1->set_right_neigbour(v2, v2->nucls()[k_ - 1]);
}

void CondensedGraph::LinkVertices(Vertex* v1, Vertex* v2) {
	DEBUG("Linking vertices '" << v1->nucls().str() << "' and '"<< v2->nucls().str() <<"' and their complement")
	assert(AreLinkable(v1, v2));

	AddRightNeighbour(v1, v2);
	AddRightNeighbour(v2->complement(), v1->complement());
}

void CondensedGraph::UnLinkVertices(Vertex* v1, Vertex* v2) {
	v1->set_right_neigbour((Vertex*) NULL, v2->nucls()[k_ - 1]);
	v2->complement()->set_right_neigbour((Vertex*) NULL,
			v1->complement()->nucls()[k_ - 1]);
}

void CondensedGraph::UnLinkAll(Vertex* v) {
	vector<Vertex*> r_ns = RightNeighbours(v);
	for (vector<Vertex*>::const_iterator it = r_ns.begin(); it != r_ns.end(); ++it) {
		UnLinkVertices(v, *it);
	}
	vector<Vertex*> l_ns = LeftNeighbours(v);
	for (vector<Vertex*>::const_iterator it = l_ns.begin(); it != l_ns.end(); ++it) {
		UnLinkVertices(*it, v);
	}
}

void DFS::ProcessVertex(Vertex* v, vector<Vertex*>& stack, Handler& h) {
	if (visited_.count(v) == 0) {
		h.HandleStartVertex(v);
		visited_.insert(v);
		vector<Vertex*> desc = g_.RightNeighbours(v);
		for (size_t i = 0; i < desc.size(); ++i) {
			Vertex* descendent = desc[i];
			h.HandleEdge(v, descendent);
			stack.push_back(descendent);
		}
	}
}

void DFS::Traverse(Handler& h) {
	for (set<Vertex*>::iterator it = g_.vertices().begin(); it
			!= g_.vertices().end(); it++) {
		vector<Vertex*> stack;
		stack.push_back(*it);
		while (!stack.empty()) {
			Vertex* v = stack[stack.size() - 1];
			stack.pop_back();
			ProcessVertex(v, stack, h);
		}
	}
}

void SimpleGraphVisualizer::Visualize(const CondensedGraph& g) {
	VisHandler h(gp_);
	DFS(g).Traverse(h);
	gp_.output();
}

void ComplementGraphVisualizer::Visualize(const CondensedGraph& g) {
	ComplementVisHandler h(gp_);
	DFS(g).Traverse(h);
	gp_.output();
}

}

