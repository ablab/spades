/**
 * condensed_graph.cpp
 *
 *  Created on: Feb 21, 2011
 *      Author: sergey
 */
#include "condensedGraph.hpp"
#include "logging.hpp"

using namespace std;

namespace condensed_graph {

Vertex::Vertex(Sequence nucls) :
	nucls_(nucls) {
	fill_n(desc_, 4, (Vertex*) NULL);
	fill_n(arc_coverage_, 4, 0);
}

Vertex::Vertex(Sequence nucls, Vertex** desc) :
	nucls_(nucls) {
	memcpy(desc, desc_, 4 * sizeof(Vertex*));
	fill_n(arc_coverage_, 4, 0);
}

Vertex::~Vertex() {
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

/**
 *	renews hash for vertex and complementary
 *	todo renew not all hashes
 */
void Graph::RenewKmersHash(Vertex* v) {
	Kmer k(v->nucls());
	h_.put(k, make_pair(v, 0));
	for (size_t i = K, n = v->nucls().size(); i < n; ++i) {
		k = k << v->nucls()[i];
		h_.put(k, make_pair(v, i - K + 1));
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

Vertex* Graph::AddVertex(Sequence nucls) {
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

bool Graph::IsMergePossible(Vertex* v1, Vertex* v2) {
	return IsLast(v1) && IsFirst(v2) && v1->complement() != v2 && v1 != v2 && AreLinkable(v1, v2);
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
	DEBUG("Merging vertices with sequences: '" << v1->nucls().str() << "' and '" << v2->nucls().str() << "' and their complement")
	assert(IsMergePossible(v1, v2));

	Vertex* v = AddVertex(v1->nucls() + v2->nucls().Subseq(K - 1));
	FixIncoming(v1, v);
	FixIncoming(v2->complement(), v->complement());

	DeleteVertex(v1);
	DeleteVertex(v2);
	return v;
}

/**
 * deletes vertex and its complement
 */
void Graph::DeleteVertex(Vertex* v) {
	DEBUG("Deleting vertex with sequence '" << v->nucls().str() << "' and its complement '" << v->complement()->nucls().str() << "'")
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
			!(v1 -> nucls()));
}

void Graph::LinkVertices(Vertex* anc, Vertex* desc) {
	DEBUG("Linking vertices with sequences '" << anc->nucls().str() << "' and '"<< desc->nucls().str() <<"'")
	assert(AreLinkable(anc, desc));

	anc->AddDesc(desc);
	component_roots_.erase(desc);
	desc->complement()->AddDesc(anc->complement());
	component_roots_.erase(anc->complement());
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
}

//pos exclusive! (goes into second vertex)
//deletes vertex if actual split happens
//returns first of the new vertices
Vertex* Graph::SplitVertex(Vertex* v, size_t pos) {
	DEBUG("Splitting vertex with sequence '" << v->nucls().str() <<"' of size " << v->size() << " at position "<< pos);
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

	DeleteVertex(v);

	return v1;
}

pair<Vertex*, int> Graph::GetPosMaybeMissing(Kmer k) {
	if (!h_.contains(k)) {
		AddVertex(Sequence(k));
	}
	return h_.get(k);
}

void Graph::ThreadRead(Read r) {
	Kmer k(r);
	DEBUG("Threading k-mer: " + k.str())
	for (size_t i = K; i < N; ++i) {
		pair<Vertex*, int> prev_pos = GetPosMaybeMissing(k);
		k = k << r[i];
		DEBUG("Threading k-mer: " + k.str())
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
					SplitVertex(v1, prev_offset + K),
					SplitVertex(v2->complement(), v2->size() - curr_offset)->complement());
		}
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

