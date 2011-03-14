#include <iostream>
#include <fstream>
#include <assert.h>
#include "logging.hpp"
#include "abruijngraph.hpp"
#include "graphVisualizer.hpp"

LOGGER("a.graph");

using namespace abruijn;

Vertex* Graph::createVertex(const Sequence* kmer) {
	LOG_ASSERT(!hasVertex(kmer), "already contains " << kmer->str());
	Vertex* v = new Vertex(kmer);
	v->complement_ = new Vertex(new Sequence(!(*kmer)));
	v->complement_->complement_ = v;
	seqVertice[*kmer] = v;
	vertices.insert(v);
	vertices.insert(v->complement_);
	return v;
}

void Graph::addEdge(Vertex* from, Vertex* to, int len) {
	from->addEdge(to, len);
	to->complement_->addEdge(from->complement_, len);
//	Edge* e = new Edge(&to, len);
//	from.addEdge(e);
//	e = new Edge(from.complement_, len);
//	to.complement_->addEdge(e);
}

void Graph::removeVertex(Vertex* v) {
	removeVertex_single(v);
	removeVertex_single(v->complement_);
}
void Graph::removeVertex_single(Vertex* v) {
	for (Vertex::Edges::iterator it = v->edges_.begin(); it != v->edges_.end(); ++it) {
		it->first->complement_->edges_.erase(v->complement_);
	}
	vertices.erase(v);
}

bool Graph::hasVertex(const Sequence* kmer) {
	return seqVertice.count(*kmer);
}

Vertex* Graph::getVertex(const Sequence* kmer) {
	LOG_ASSERT(hasVertex(kmer), "No such vertex: " << kmer->str());
	Vertex* v = seqVertice[*kmer];
	if (*kmer == *v->kmer_) {
		return v;
	}
	LOG_ASSERT(*kmer == *v->complement_->kmer_, "Bug in reverse-complimentary");
	return v->complement_;
}

static int k = 0;

Vertex* Graph::condense(Vertex* v) {
	Vertex* u = v->edges_.begin()->first;
	Edge e = v->edges_.begin()->second;
	assert(u->complement_->edges_.begin()->first == v->complement_);
	Edge f = u->complement_->edges_.begin()->second;
	DEBUG(e.toString() << " " << f.toString());
	if ((e.lengths_.size() != 1) || (e.lengths_.size() != 1)) {
		ERROR("Unsure what to do");
		return NULL;
	}
	int len = e.lengths_.begin()->first;
	assert(len == f.lengths_.begin()->first);
	SequenceBuilder sb;
	sb.append(*(v->kmer_));
	if (len > v->kmer_->size()) {
		DEBUG("Should've been filled correctly.");
		for (int i = v->kmer_->size(); i < len; i++) {
			sb.append('A');
		}
		sb.append(*(u->kmer_));
	} else {
		sb.append(u->kmer_->Subseq(v->kmer_->size() - len));
	}
	Vertex* vu = createVertex(new Sequence(sb.BuildSequence()));
	for (Vertex::Edges::iterator it = u->edges_.begin(); it != u->edges_.end(); ++it) {
		Vertex* w = it->first;
		for (map<int, int>::iterator it2 = it->second.lengths_.begin(); it2 != it->second.lengths_.end(); ++it2) {
			vu->edges_[w].lengths_[len + it2->first] += it2->second;
			w->complement_->edges_[vu->complement_].lengths_[len + it2->first] += it2->second;
		}
	}
	for (Vertex::Edges::iterator it = v->complement_->edges_.begin(); it != v->complement_->edges_.end(); ++it) {
		Vertex* w = it->first;
		for (map<int, int>::iterator it2 = it->second.lengths_.begin(); it2 != it->second.lengths_.end(); ++it2) {
			vu->complement_->edges_[w].lengths_[len + it2->first] += it2->second;
			w->complement_->edges_[vu].lengths_[len + it2->first] += it2->second;
		}
	}
	k++;
	removeVertex(v);
	removeVertex(u);
	return vu;
}

const int B = 3;
string shorten(const Sequence* s) {
	return s->Subseq(0, B).str() + "_" + itoa(s->size() - 2 * B) + "_"+ s->Subseq(s->size() - B).str();
}

void Graph::output(std::ofstream &out) {
	gvis::GraphPrinter<string> printer("z", out);
	for (Vertices::iterator v = vertices.begin(); v != vertices.end(); ++v) {
		printer.addVertex((*v)->kmer_->str(), shorten((*v)->kmer_));
//		for (Vertex::Edges::iterator e = (*v)->edges_.begin(); e != (*v)->edges_.end(); ++e) {
//			printer.addEdge((*v)->kmer_->str(), e->to_->kmer_->str(), "oppa"); // TODO len
//		}
		for (Vertex::Edges::iterator it = (*v)->edges_.begin(); it != (*v)->edges_.end(); ++it) {
			stringstream ss;
			printer.addEdge((*v)->kmer_->str(), it->first->kmer_->str(), it->second.toString());
		}
	}
	printer.output();
	out.close();
}

void Graph::output(string filename) {
	ofstream out(filename.c_str(), ios::out);
	output(out);
}

