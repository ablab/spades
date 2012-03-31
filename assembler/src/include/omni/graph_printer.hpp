#pragma once

#include "graph_print_utils.hpp"

namespace omnigraph {
//todo ask Valera if it is ok
using namespace gvis;

template<typename VertexId>
class GraphPrinter {
protected:
	const string name_;
	ostream& out_;
public:
	GraphPrinter(const string &name, ostream &out = cout) :
		name_(name), out_(out) {
	}

	virtual ~GraphPrinter() {
	}

	virtual void open() = 0;

	virtual void close() = 0;

	virtual void AddVertex(VertexId vertexId, const string &label = " ",
			const string fillColor = "white") = 0;

	virtual void AddEdge(VertexId fromId, VertexId toId,
			const string &label = " ", const string color = "black", const int length = 0) = 0;

};

template<typename VertexId>
class DotGraphPrinter: public GraphPrinter<VertexId> {
private:
	typedef GraphPrinter<VertexId> super;
public:
	DotGraphPrinter(const string &name, ostream &out = cout) :
		super(name, out) {

	}

	virtual void open() {
		startSimpleGraphRecord(super::out_, super::name_);
	}

	virtual void AddVertex(VertexId vertexId, const string &label,
			const string fillColor = "white") {
		Vertex<VertexId> v(vertexId, label, fillColor);
		recordVertex<VertexId> (super::out_, v);
	}

	virtual void AddEdge(VertexId fromId, VertexId toId, const string &label,
			const string color = "black", const int length = 0) {
		Edge<VertexId> e(fromId, toId, label, color, length);
		recordSimpleEdge<VertexId> (super::out_, e);
	}

	virtual void close() {
		endGraphRecord(super::out_);
	}

	virtual ~DotGraphPrinter() {

	}

};

//todo now it is even worse than ever!!!
template<class Graph, typename VertexId = typename Graph::VertexId>
class DotPairedGraphPrinter: public GraphPrinter<VertexId> {
private:
	typedef GraphPrinter<VertexId> super;
	const Graph& g_;
	const GraphLabeler<Graph> &labeler_;
	PairedGraphPrinter<VertexId> paired_printer_;
public:
	DotPairedGraphPrinter(const Graph& g, const string &name, const GraphLabeler<Graph> &labeler,
			ostream &out = cout) :
		super(name, out), g_(g), labeler_(labeler), paired_printer_(name, out) {
	}

	virtual ~DotPairedGraphPrinter() {

	}

	virtual void open() {
		paired_printer_.open();
	}

	virtual void close() {
		paired_printer_.close();
	}

	virtual void AddVertex(VertexId v, const string &label,
			const string fillColor) {
		string label_conjugate = labeler_.label(g_.conjugate(v));
		paired_printer_.AddVertex(v, label, g_.conjugate(v), label_conjugate, fillColor);
	}

	virtual void AddEdge(VertexId v1, VertexId v2, const string &label,
			const string color, const int length = 0) {
		paired_printer_.AddEdge(make_pair(v1, g_.conjugate(v1)),
				make_pair(v2, g_.conjugate(v2)), label, color, length);
	}
};

}
