/**
 * condensed_graph.cpp
 *
 *  Created on: Feb 21, 2011
 *      Author: sergey
 */
#include "condensedGraph.hpp"
#include "logging.hpp"
#include "graphVisualizer.hpp"
#include <set>
#include <iostream>
#include <tr1/unordered_map>

using namespace std;

namespace condensed_graph {

const set<Vertex*>& Graph::component_roots() const {
	return component_roots_;
}

vector<Vertex*> Graph::Anc(const Vertex* v) const {
	vector<Vertex*> ans;
	Vertex* const * compl_desc = v->complement()->desc();
	for (int i = 3; i >= 0; --i) {
		if (compl_desc[i] != NULL) {
			ans.push_back(compl_desc[i]->complement());
		}
	}
	return ans;
}

vector<Vertex*> Graph::Desc(const Vertex* v) const {
	vector<Vertex*> ans;
	Vertex* const * desc = v->desc();
	for (int i = 0; i < 4; ++i) {
		if (desc[i] != NULL) {
			ans.push_back(desc[i]);
		}
	}
	return ans;
}

/*
 bool Graph::AddIfRoot(Vertex* v) {
 bool f = false;
 if (Anc(v).empty()) {
 component_roots_.insert(v);
 f = true;
 }
 if (Desc(v).empty()) {
 component_roots_.insert(v->complement());
 f = true;
 }
 return f;
 }
 */

/**
 * adds vertex and its complement
 *///std::string operator+(const std::string& s, int i) {
//	std::stringstream out;
//	out << s << i;
//	return out.str();
//}

Vertex* Graph::AddVertex(const Sequence &nucls) {
	DEBUG("Adding vertex for sequence '" << nucls.str() << "' and its complement '" << (!nucls).str() << "'")
	Vertex* v1 = new Vertex(nucls);
	Vertex* v2 = new Vertex(!nucls);
	v1->set_complement(v2);
	v2->set_complement(v1);
	component_roots_.insert(v1);
	component_roots_.insert(v2);
	//	DEBUG("Renewing hash for k-mers of sequence " << v->nucls().str() << " and its complement")
	RenewKmersHash(v1);
	RenewKmersHash(v2);
	return v1;
}

bool Graph::IsMergePossible(Vertex* v1, Vertex* v2) const {
	return IsLast(v1) && IsFirst(v2) && v1->complement() != v2 && v1 != v2
			&& AreLinkable(v1, v2);
}

bool Graph::CanBeDeleted(Vertex* v) const {
	vector<Vertex*> anc = Anc(v);
	for (size_t i = 0; i < anc.size(); ++i) {
		Vertex* ancestor = anc[i];
		if (ancestor != v && ancestor != v->complement()) {
			for (size_t j = 0; j < 4; ++j) {
				if (ancestor->desc()[j] == v) {
					return false;
				}
			}
		}
	}
	return true;
}

/**
 * deals with incoming links and their complement only!!!
 */
void Graph::FixIncomingOnSplit(Vertex* v, Vertex* v1, Vertex* v2) {
	vector<Vertex*> anc = Anc(v);
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

//pos exclusive! (goes into second vertex)
//deletes vertex if actual split happens
//returns first of the new vertices
Vertex* Graph::SplitVertex(Vertex* v, size_t pos) {
	DEBUG("Splitting vertex '" << v->nucls().str() <<"' of size " << v->size() << " at position "<< pos);
	assert(pos <= v->size());

	if (pos == v->size()) {
		return v;
	};

	Sequence nucls = v->nucls();

	Vertex* v1 = AddVertex(nucls.Subseq(0, pos));
	Vertex* v2 = AddVertex(nucls.Subseq(pos - (K - 1), nucls.size())); // nucls.size() can be omitted here

	LinkVertices(v1, v2);

	FixIncomingOnSplit(v, v1, v2);

	FixIncomingOnSplit(v->complement(), v2->complement(), v1->complement());

	DeleteVertex(v);

	return v1;
}

void Graph::FixIncomingOnMerge(Vertex* v1, Vertex* v2, Vertex* v) {
	vector<Vertex*> anc = Anc(v1);
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

Vertex* Graph::Merge(Vertex* v1, Vertex* v2) {
	DEBUG("Merging vertices '" << v1->nucls().str() << "' and '" << v2->nucls().str() << "' and their complement")
	assert(IsMergePossible(v1, v2));

	Vertex* v = AddVertex(v1->nucls() + v2->nucls().Subseq(K - 1));
	FixIncomingOnMerge(v1, v2, v);
	FixIncomingOnMerge(v2->complement(), v1->complement(), v->complement());

	DeleteVertex(v1);
	DeleteVertex(v2);
	return v;
}

/**
 * deletes vertex and its complement
 */
void Graph::DeleteVertex(Vertex* v) {
	DEBUG("Deleting vertex '" << v->nucls().str() << "' and its complement '" << v->complement()->nucls().str() << "'")
	assert(CanBeDeleted(v));
	assert(CanBeDeleted(v->complement()));

	Vertex* complement = v->complement();
	component_roots_.erase(v);
	component_roots_.erase(complement);

	delete v;
	delete complement;
}

bool Graph::AreLinkable(Vertex* v1, Vertex* v2) const {
	return KMinusOneMer(v2->nucls()) == KMinusOneMer(
			v1->nucls(),
			v1->size() - (K - 1));// && !v1->deleted && !v2 -> deleted;
	//was: return KMinusOneMer(v2 -> nucls()) == !KMinusOneMer(!v1 -> nucls());// && !v1->deleted && !v2 -> deleted;
}

void Graph::LinkVertices(Vertex* anc, Vertex* desc) {
	DEBUG("Linking vertices '" << anc->nucls().str() << "' and '"<< desc->nucls().str() <<"' and their complement")
	assert(AreLinkable(anc, desc));

	anc->AddDesc(desc);
	//component_roots_.erase(desc);
	desc->complement()->AddDesc(anc->complement());
	//component_roots_.erase(anc->complement());
}

void Graph::ThreadRead(const Read &r) {
	Kmer k(r);
	DEBUG("Threading k-mer: " + k.str())
	for (size_t i = K; i < N; ++i) {
		pair<Vertex*, int> prev_pos = GetPosMaybeMissing(k);
		Kmer old_k = k;
		k = k << r[i];
		DEBUG("Threading k-mer: " + k.str())
		pair<Vertex*, int> curr_pos = GetPosMaybeMissing(k);

		Vertex* prev_v = prev_pos.first;
		Vertex* curr_v = curr_pos.first;
		size_t prev_offset = prev_pos.second;
		size_t curr_offset = curr_pos.second;

		if (IsLastKmer(prev_v, prev_offset) && IsFirstKmer(curr_v, curr_offset)
				&& IsMergePossible(prev_v, curr_v)) {
			Merge(prev_v, curr_v);
		} else if (prev_v == curr_v && prev_offset + 1 == curr_offset) {
			//todo check links here to optimize???
			//do nothing
		} else {
			SplitVertex(prev_v, prev_offset + K);
			//need if k-mers were on same or complementary vertices
			curr_pos = GetPosition(k);
			Vertex* curr_v = curr_pos.first;
			size_t curr_offset = curr_pos.second;
			Vertex* v2 = SplitVertex(curr_v->complement(),
					curr_v->size() - curr_offset)->complement();
			Vertex* v1 = GetPosition(old_k).first;
			LinkVertices(v1, v2);
		}
	}
}

Traversal::Traversal(const Graph& g) :
	g_(g) {
}

DFS::DFS(const Graph& g) :
	Traversal(g) {
}

void DFS::Traverse(Handler& h) {
	for (set<Vertex*>::iterator it = g_.component_roots().begin(); it != g_.component_roots().end(); it++) {
		vector<Vertex*> stack;
		stack.push_back(*it);
		while (!stack.empty()) {
			go(stack[stack.size() - 1], stack, h);
			stack.pop_back();
		}
	}
}

void DFS::go(Vertex* v, vector<Vertex*>& stack, Handler& h) {
	if (visited_.count(v) == 0) {
		h.HandleStartVertex(v);
		visited_.insert(v);
		vector<Vertex*> desc = g_.Desc(v);
		for (size_t i = 0; i < desc.size(); ++i) {
			Vertex* descendent = desc[i];
			h.HandleEdge(v, descendent);
			stack.push_back(descendent);
		}
	}
}

class CountHandler: public Traversal::Handler {
	tr1::unordered_map<const Vertex*, size_t>& map_;
	size_t count_;
public:

	CountHandler(tr1::unordered_map<const Vertex*, size_t>& map) :
		map_(map), count_(0) {
	}

	virtual void HandleStartVertex(const Vertex* v) {
		map_.insert(make_pair(v, count_++));
	}
};

class VisHandler: public Traversal::Handler {
	gvis::IGraphPrinter<Vertex*>& pr_;
public:

	VisHandler(gvis::IGraphPrinter<Vertex*>& pr) :
		pr_(pr) {
	}

	virtual void HandleStartVertex(Vertex* v) {
		stringstream ss;
		ss << v->nucls().size();
		pr_.addVertex(v, ss.str());
	}

	virtual void HandleEdge(Vertex* v1, Vertex* v2) {
		pr_.addEdge(v1, v2, "");
	}

};

void SimpleGraphVisualizer::Visualize(const Graph& g) {
	VisHandler h(gp_);
	DFS(g).Traverse(h);
	gp_.output();
}

}

