#include <iostream>
#include <fstream>
#include <assert.h>
#include "logging.hpp"
#include "abruijngraph.hpp"
#include "graphVisualizer.hpp"
#include "simple_tools.hpp"
#include <ext/hash_map>

namespace abruijn {

ostream& operator<< (ostream& os, const Edge& e) {
	return e.output(os);
}

ostream& operator<< (ostream& os, const Vertex& v) {
	return v.output(os);
}

Vertex* Graph::createVertex(const Sequence& kmer) {
	LOG_ASSERT(!hasVertex(kmer), "already contains " << kmer.str());
	Vertex* v = new Vertex(kmer, true);
	seqVertice[kmer] = v;
	seqVertice[!kmer] = v->complement();
	vertices.insert(v);
	vertices.insert(v->complement());
	return v;
}

//Vertex* Graph::createVertex(const Sequence& kmer, size_t size) {
//	Vertex* v = createVertex(kmer);
//	v->setSize(size);
//	v->complement()->setSize(size);
//	return v;
//}


void Graph::addEdge(Vertex* from, Vertex* to, const Sequence& seq) {
	from->addEdge(to, seq);
	to->complement()->addEdge(from->complement(), !seq);
}

void Graph::removeVertex(Vertex* v) {
	removeVertex_single(v);
	removeVertex_single(v->complement());
}
void Graph::removeVertex_single(Vertex* v) {
	for (Edges::iterator it = v->edges().begin(); it != v->edges().end(); ++it) {
		it->first->complement()->edges().erase(v->complement());
	}
	removedVertices.insert(v);
}

bool Graph::hasVertex(const Sequence& kmer) const {
	return seqVertice.count(kmer);
}

/**
 * @return Find or create a vertex corresponding to this kmer
 */
Vertex* Graph::getVertex(const Sequence& kmer) {
	SeqVertice::iterator v = seqVertice.find(kmer);
	if (v == seqVertice.end()) {
		return createVertex(kmer);
	}
	if (v->second->data() == kmer) {
		return v->second;
	}
	assert(v->second->complement()->data() == kmer);
	return v->second->complement();
}

void Graph::Condense_single(Vertex* v) {
	assert(v->degree() == 1 && v->complement()->degree() == 1);
	Vertex* u = v->edges().begin()->first;
	Vertex* w = v->complement()->edges().begin()->first->complement();
	assert(u != NULL);
	assert(w != NULL);
//	DEBUG(w->edges().count(u));
//	assert(w->edges().count(u) == 0);
	Edge e1 = w->edges()[v];
	Edge e2 = v->edges()[w];
	Edge e3 = e1.concat(e2);

//	(w->edges())[u] = Edge();
//	(w->edges())[u] = (w->edges())[u];
//	cout << ((w->edges())[u] = Edge()) <<endl;
//	w->edges()[u] = v->edges()[w].add(w->edges()[u]);
	w->edges()[u] += e3;
}

bool Graph::Condense(Vertex* v) {
	if (v->degree() != 1 || v->complement()->degree() != 1) {
		return false;
	}
//	if (v->backward()->edges().count(v->forward()) != 0) {
//		return false;
//	}
//	if (v->forward()->complement()->edges().count(v->backward()->complement()) != 0) {
//		return false;
//	}
	assert(v->backward()->edges().count(v->forward()) == 0);
	assert(v->forward()->complement()->edges().count(v->backward()->complement()) == 0);
//	v->backward()->edges()[v->forward()].touch();
//	v->forward()->complement()->edges()[v->backward()->complement()].touch();
//	DEBUG("touched");
//	Condense_single(v);
//	Condense_single(v->complement());
//	removeVertex(v);
	return false;
//	Vertex* w = v->complement()->edges_.begin()->first->complement();
//	Vertex* u = v->edges_.begin()->first;
//	for (;;) {
//
//	}

//	Edge e = v->edges_.begin()->second;
//	assert(u->complement()->edges_.begin()->first == v->complement());
//	Edge f = u->complement()->edges_.begin()->second;
//	TRACE(e << " " << f);
//	if ((e.lengths_.size() != 1) || (e.lengths_.size() != 1)) {
//		ERROR("Unsure what to do");
//		return NULL;
//	}
//	size_t len = e.lengths_.begin()->first;
//	assert(len == f.lengths_.begin()->first);
//	Vertex* vu = createVertex(v->concat(u), len);
//	for (Vertex::Edges::iterator it = u->edges_.begin(); it != u->edges_.end(); ++it) {
//		Vertex* w = it->first;
//		for (map<size_t, int>::iterator it2 = it->second.lengths_.begin(); it2 != it->second.lengths_.end(); ++it2) {
//			size_t l = len + it2->first - u->size();
//			vu->edges_[w].lengths_[l] += it2->second;
//			w->complement()->edges_[vu->complement()].lengths_[l] += it2->second;
//		}
//	}
//	for (Vertex::Edges::iterator it = v->complement()->edges_.begin(); it != v->complement()->edges_.end(); ++it) {
//		Vertex* w = it->first;
//		for (map<size_t, int>::iterator it2 = it->second.lengths_.begin(); it2 != it->second.lengths_.end(); ++it2) {
//			size_t l = len + it2->first - v->size();
//			vu->complement()->edges_[w].lengths_[l] += it2->second;
//			w->complement()->edges_[vu].lengths_[l] += it2->second;
//		}
//	}
//	k++;
//	removeVertex(v);
//	removeVertex(u);
//	return vu;
}

void Graph::Condense() {
	DEBUG(vertices.size() << " vertices");
	int condensations = 0;
	for (;;) {
		bool change = false;
		for (iterator v = begin(); v != end(); ++v) {
			if (Condense(*v)) {
				condensations++;
				change = true;
			}
		}
		cleanup();
		INFO(condensations << " condensations");
		if (!change) {
			break;
		}
	}
	DEBUG(vertices.size() << " vertices");
}

void Graph::cleanup() {
	for (Vertices::iterator v = removedVertices.begin(); v != removedVertices.end(); ++v) {
		vertices.erase(*v);
		delete *v;
	}
	removedVertices.clear();
}

void Graph::stats() {
	size_t tips = 0;
	for (iterator v = begin(); v != end(); ++v) {
		if (v->degree() == 0) {
			tips++;
		}
	}
	INFO(vertices.size() << " vertices");
	INFO(tips << " tips");
}

void Graph::output(std::ofstream &out, bool paired) {
	std::string name = "A_Bruijn_Graph";
	if (paired) {
		gvis::PairedGraphPrinter<Vertex*> printer(name, out);
		for (Vertices::iterator v = vertices.begin(); v != vertices.end(); ++v) {
			if ((*v)->data() < (*v)->complement()->data()) {
				printer.addVertex(*v, ToString(**v), (*v)->complement(), ToString(*((*v)->complement())));
			}
			for (Edges::iterator it = (*v)->edges().begin(); it != (*v)->edges().end(); ++it) {
				printer.addEdge(make_pair(*v, (*v)->complement()), make_pair(it->first, it->first->complement()), ToString(it->second));
			}
		}
		printer.output();
	} else {
		gvis::GraphPrinter<Vertex*> printer(name, out);
		for (Vertices::iterator v = vertices.begin(); v != vertices.end(); ++v) {
			printer.addVertex(*v, ToString(**v));
			for (Edges::iterator it = (*v)->edges().begin(); it != (*v)->edges().end(); ++it) {
				printer.addEdge(*v, it->first, ToString(it->second));
			}
		}
	}
}

void Graph::output(string filename, bool paired) {
	ofstream out(filename.c_str(), ios::out);
	output(out, paired);
	out.close();
}

}
