/**
 * condensed_graph.cpp
 *
 *  Created on: Feb 21, 2011
 *      Author: sergey
 */
#include "condensed_graph.hpp"
#include "logging.hpp"

using namespace std;

LOGGER("debruijn.condensed_graph")

namespace condensed_graph {

Vertex::Vertex(Sequence nucls) :
	nucls_(nucls) {
	fill_n(desc_, 4, (Vertex*) NULL);
}

//todo talk with Kolya about problem with Sequence copying!!!
Vertex::Vertex(Sequence nucls, Vertex** desc) :
	nucls_(nucls) {
	memcpy(desc, desc_, 4 * sizeof(Vertex*));
	//		_arcs = arcs;
}

Vertex::~Vertex() {
	delete[] desc_;
	delete[] arc_coverage_;
}

int Vertex::DescCount() {
	int c = 0;
	for (int i = 0; i < 4; ++i)
		if (desc_[i] != NULL)
			c++;
	return c;
}

bool Vertex::IsDeadend() {
	return DescCount() == 0;
}

Vertex** Vertex::desc() {
	return desc_;
}

Vertex* Vertex::desc(char nucl) {
	return desc_[(int) nucl];
}

size_t Vertex::size() {
	return nucls_.size();
}

Sequence Vertex::nucls() {
	return nucls_;
}

void Vertex::AddDesc(Vertex* v) {
	//check if k-1 mers differ
	int k_th_nucl = v->nucls()[K - 1];
	desc_[k_th_nucl] = v;
}

Vertex* Vertex::complement() {
	return complement_;
}

void Vertex::set_complement(Vertex* complement) {
	complement_ = complement;
}

int Vertex::coverage() {
	return coverage_;
}

void Vertex::set_coverage(int coverage) {
	coverage_ = coverage;
}

void Graph::RenewHashForSingleVertexKmers(Vertex* v) {
	Kmer k(v->nucls());
	h_.put(k, make_pair(v, 0));
	for (size_t i = K + 1, n = v->nucls().size(); i < n; ++i) {
		k = k << v->nucls()[i];
		h_.put(k, make_pair(v, i - K));
	}
}

const set<Vertex*>& Graph::component_roots() {
	return component_roots_;
}

vector<Vertex*> Graph::Anc(Vertex* v) {
	vector<Vertex*> ans;
	Vertex** compl_desc = v->complement()->desc();
	for (int i = 3; i >= 0; --i) {
		if (compl_desc[i] != NULL) {
			ans.push_back(compl_desc[i]->complement());
		}
	}
	return ans;
}

vector<Vertex*> Graph::Desc(Vertex* v) {
	vector<Vertex*> ans;
	Vertex** desc = v->desc();
	for (int i = 0; i < 4; ++i) {
		if (desc[i] != NULL) {
			ans.push_back(desc[i]);
		}
	}
	return ans;
}

bool Graph::IsLast(Vertex* v) {
	return v->IsDeadend();
}

bool Graph::IsFirst(Vertex* v) {
	return IsLast(v->complement());
}

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

/**
 * adds vertex and its complement
 */
Vertex* Graph::AddVertex(Sequence nucls) {
	Vertex* v1 = new Vertex(nucls);
	Vertex* v2 = new Vertex(!nucls);
	v1->set_complement(v2);
	v2->set_complement(v1);
	return v1;
}

bool Graph::IsMergePossible(Vertex* v1, Vertex* v2) {
	return IsLast(v1) && IsFirst(v2) && v1->complement() != v2 && v1 != v2;
}

bool Graph::CanBeDeleted(Vertex* v) {
	vector<Vertex*> anc = Anc(v);
	for (size_t i = 0; i < anc.size(); ++i) {
		for (size_t j = 0; j < 4; ++j) {
			if (anc[i]->desc()[j] == v) {
				return false;
			}
		}
	}
	return true;
}

Vertex* Graph::Merge(Vertex* v1, Vertex* v2) {
	assert(IsMergePossible(v1, v2));

	Vertex* v = AddVertex(v1->nucls() + v2->nucls());
	FixIncoming(v1, v);
	FixIncoming(v2->complement(), v->complement());
	RenewHashForVertexKmers(v);
	return v;
}

/**
 * deletes vertex and its complement
 */
void Graph::DeleteVertex(Vertex* v) {
	assert(CanBeDeleted(v));
	assert(CanBeDeleted(v->complement()));

	Vertex* complement = v->complement();
	component_roots_.erase(v);
	component_roots_.erase(complement);
	delete v;
	delete complement;
}

bool Graph::AreLinkable(Vertex* v1, Vertex* v2) {
	return KMinusOneMer(v2 -> nucls()) == !KMinusOneMer(
			v1 -> complement() -> nucls());
}

void Graph::LinkVertices(Vertex* anc, Vertex* desc) {
	assert(AreLinkable(anc, desc));

	anc->AddDesc(desc);
	component_roots_.erase(anc);
	desc->complement()->AddDesc(anc->complement());
	component_roots_.erase(desc->complement());
}

/**
 * deals with incoming links and their complement only!!!
 */
void Graph::FixIncoming(Vertex* v, Vertex* new_v) {
	assert(Kmer(v->nucls()) == Kmer(v->nucls()));

	vector<Vertex*> anc = Anc(v);
	for (size_t i = 0; i < anc.size(); ++i) {
		LinkVertices(anc[i], new_v);
	}
	AddIfRoot(new_v);
}

/**
 *	renews hash for vertex and complementary
 *	todo renew not all hashes
 */
void Graph::RenewHashForVertexKmers(Vertex* v) {
	RenewHashForSingleVertexKmers(v);
	RenewHashForSingleVertexKmers(v->complement());
}

//pos exclusive! (goes into second vertex)
//deletes vertex if actual split happens
Vertex* Graph::SplitVertex(Vertex* v, size_t pos, bool return_second) {
	assert(pos <= v->size());

	if (pos == v->size()) {
		return v;
	};

	Sequence nucls = v->nucls();

	Vertex* v1 = AddVertex(nucls.Subseq(0, pos));
	Vertex* v2 = AddVertex(nucls.Subseq(pos - (K - 1), nucls.size()));

	LinkVertices(v1, v2);

	FixIncoming(v, v1);

	FixIncoming(v->complement(), v2->complement());

	RenewHashForVertexKmers(v1);

	RenewHashForVertexKmers(v2);

	DeleteVertex(v);

	return return_second ? v2 : v1;
}

pair<Vertex*, int> Graph::GetPosMaybeMissing(Kmer k) {
	if (!h_.contains(k)) {
		AddVertex(Sequence(k));
	}
	return h_.get(k);
}

void Graph::ThreadRead(Read r) {
	Kmer k(r);
	pair<Vertex*, int> prev_pos = GetPosMaybeMissing(k);
	for (size_t i = K + 1; i < N; ++i) {
		k = k << r[i];
		pair<Vertex*, int> curr_pos = GetPosMaybeMissing(k);

		Vertex* v1 = prev_pos.first;
		Vertex* v2 = curr_pos.first;
		size_t prev_offset = prev_pos.second;
		size_t curr_offset = curr_pos.second;

		if (IsLastKmer(v1, prev_offset) && IsFirstKmer(v2, curr_offset)
				&& IsMergePossible(v1, v2)) {
			Merge(v1, v2);
		} else if (v1 == v2 && prev_offset + 1 == curr_offset) {
			//do nothing
		} else {
			LinkVertices(
					SplitVertex(v1, prev_offset + K, true),
					SplitVertex(v2->complement(), N - curr_offset, false)->complement());
		}
		prev_pos = curr_pos;
	}
}

bool Graph::IsLastKmer(Vertex* v, size_t pos) {
	return pos + K == v->size();
}

bool Graph::IsFirstKmer(Vertex* v, size_t pos) {
	return pos == 0;
}



/*VertexPool::VertexPool(_v_idx size) : _size(size) {
	_max_idx = 0;
	_free = new bool[size];
	fill_n(_free, size, true);
	_vertices = new Vertex[size];
}

VertexPool::~VertexPool() {
	delete [] _free;
	delete [] _vertices;
}

const Vertex& operator[](_v_idx index) const;
bool IsFree(_v_idx index) const;
void Free(_v_idx index);
_v_idx AllocatePair();
}*/
}

