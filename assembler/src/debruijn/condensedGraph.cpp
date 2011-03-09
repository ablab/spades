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

Vertex::Vertex(const Sequence &nucls) :
	nucls_(nucls) {
	fill_n(desc_, 4, (Vertex*) NULL);
	fill_n(arc_coverage_, 4, 0);
//	deleted = false;
}

Vertex::Vertex(const Sequence &nucls, Vertex** desc) :
	nucls_(nucls) {
	memcpy(desc, desc_, 4 * sizeof(Vertex*));
	fill_n(arc_coverage_, 4, 0);
//	deleted = false;
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
	h_.put(k, v, 0);
	for (size_t i = K, n = v->nucls().size(); i < n; ++i) {
		k = k << v->nucls()[i];
		h_.put(k, v, i - K + 1);
	}
}

const set<Vertex*>& Graph::component_roots() const {
	return component_roots_;
}

vector<Vertex*> Graph::Anc(Vertex* v) const {
	vector<Vertex*> ans;
	Vertex** compl_desc = v->complement()->desc();
	for (int i = 3; i >= 0; --i) {
		if (compl_desc[i] != NULL) {
			ans.push_back(compl_desc[i]->complement());
		}
	}
	return ans;
}

vector<Vertex*> Graph::Desc(Vertex* v) const {
	vector<Vertex*> ans;
	Vertex** desc = v->desc();
	for (int i = 0; i < 4; ++i) {
		if (desc[i] != NULL) {
			ans.push_back(desc[i]);
		}
	}
	return ans;
}

bool Graph::IsLast(Vertex* v) const {
	return v->IsDeadend();
}

bool Graph::IsFirst(Vertex* v) const {
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
	return true;// && !v->deleted;
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
	Vertex* v2 = AddVertex(nucls.Subseq(pos - (K - 1), nucls.size()));

	LinkVertices(v1, v2);

	FixIncomingOnSplit(v, v1, v2);

	FixIncomingOnSplit(v->complement(), v2->complement(), v1 -> complement());

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

//	v->deleted = true;
//	complement->deleted = true;
	delete v;
	delete complement;
}

bool Graph::AreLinkable(Vertex* v1, Vertex* v2) const {
	return KMinusOneMer(v2 -> nucls()) == !KMinusOneMer(!(v1 -> nucls()));// && !v1->deleted && !v2 -> deleted;
}

void Graph::LinkVertices(Vertex* anc, Vertex* desc) {
	DEBUG("Linking vertices '" << anc->nucls().str() << "' and '"<< desc->nucls().str() <<"' and their complement")
	assert(AreLinkable(anc, desc));

	anc->AddDesc(desc);
	component_roots_.erase(desc);
	desc->complement()->AddDesc(anc->complement());
	component_roots_.erase(anc->complement());
}

pair<Vertex*, int> Graph::GetPosMaybeMissing(Kmer k) {
	if (!h_.contains(k)) {
		AddVertex(Sequence(k));
	}
	return h_.get(k);
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
			//need if k-mers were on same complementary vertices
			curr_pos = GetPosMaybeMissing(k);
			Vertex* curr_v = curr_pos.first;
			size_t curr_offset = curr_pos.second;
			Vertex* v2 = SplitVertex(curr_v->complement(), curr_v->size() - curr_offset)->complement();
			Vertex* v1 = GetPosMaybeMissing(old_k).first;
			LinkVertices(v1, v2);
		}
	}
}

bool Graph::IsLastKmer(Vertex* v, size_t pos) const {
	return pos + K == v->size();
}

bool Graph::IsFirstKmer(Vertex* v, size_t pos) const {
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

