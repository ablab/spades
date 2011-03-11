#include <iostream>
#include <fstream>
#include "logging.hpp"
#include "abruijngraph.hpp"
#include "graphVisualizer.hpp"

LOGGER("a.graph");

using namespace abruijn;

void Graph::addVertex(const Sequence* kmer) {
	LOG_ASSERT(!hasVertex(kmer), "already contains " << kmer->str());
	Vertex* v = new Vertex(kmer);
	v->complement_ = new Vertex(new Sequence(!(*kmer)));
	v->complement_->complement_ = v;
	seqVertice[*kmer] = v;
	vertexArray.push_back(v);
	vertexArray.push_back(v->complement_);
}

void Graph::addEdge(Vertex &from, Vertex &to, int len) {
	TRACE(&from);
	Edge* e = new Edge(&to, len);
	from.addEdge(e);
	TRACE(from.edges_.size());
	e = new Edge(from.complement_, len);
	to.complement_->addEdge(e);
}

bool Graph::hasVertex(const Sequence* kmer) {
	return seqVertice.count(*kmer);
}

abruijn::Vertex* Graph::getVertex(const Sequence* kmer) {
	LOG_ASSERT(hasVertex(kmer), "No such vertex: " << kmer->str());
	Vertex* v = seqVertice[*kmer];
	if (*kmer == *v->kmer_) {
		return v;
	}
	LOG_ASSERT(*kmer == *v->complement_->kmer_, "Bug in reverse-complimentary");
	return v->complement_;
}

void Graph::output(const char* filename) {
	std::ofstream out(filename, ios::out);
	gvis::GraphPrinter<string> printer("z", out);
	for (VertexArray::iterator v = vertexArray.begin(); v != vertexArray.end(); ++v) {
		printer.addVertex((*v)->kmer_->str(), (*v)->kmer_->str());
		for (Vertex::EdgeArray::iterator e = (*v)->edges_.begin(); e != (*v)->edges_.end(); ++e) {
			printer.addEdge((*v)->kmer_->str(), e->to_->kmer_->str(), "oppa"); // TODO len
		}
	}
	printer.output();
	out.close();
}
