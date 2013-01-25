//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "graph_print_utils.hpp"
#include "graph_labeler.hpp"
#include "graph_colorer.hpp"

namespace omnigraph {
//todo strange parameter length... is it obsollete???

//todo ask Valera if it is ok
using namespace gvis;

template<class VertexId, class EdgeId>
class GraphElementPrinter {
	const string name_;
	ostream& out_;
protected:
	const string& name() {
		return name_;
	}

	ostream& out() {
		return out_;
	}

public:
	GraphElementPrinter(const string &name, ostream &out = cout) :
		name_(name), out_(out) {
	}

	virtual ~GraphElementPrinter() {
	}

	virtual void open() = 0;

	virtual void close() = 0;

	virtual void AddVertex(VertexId v) = 0;
	virtual void AddVertices(const std::set<VertexId> &elements) = 0;

	virtual void AddEdge(EdgeId e, const int length = 0) = 0;

};

template<class Graph>
class GraphPrinter : public GraphElementPrinter<typename Graph::VertexId, typename Graph::EdgeId> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef GraphElementPrinter<VertexId, EdgeId> base;

	const Graph& g_;
	const GraphLabeler<Graph>& labeler_;
	const GraphColorer<Graph>& colorer_;
protected:

	const Graph& g() {
		return g_;
	}

	const GraphLabeler<Graph>& labeler() {
		return labeler_;
	}

	const GraphColorer<Graph>& colorer() {
		return colorer_;
	}

public:
	GraphPrinter(const Graph& g, const GraphLabeler<Graph>& labeler,
			const GraphColorer<Graph>& colorer, const string &name, ostream &out = cout) :
		base(name, out), g_(g), labeler_(labeler), colorer_(colorer) {
	}

};

template<class Graph>
class DotGraphPrinter: public GraphPrinter<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef GraphPrinter<Graph> base;
public:
	DotGraphPrinter(const Graph& g, const GraphLabeler<Graph>& labeler,
			const GraphColorer<Graph>& colorer, const string &name, ostream &out = cout) :
		base(g, labeler, colorer, name, out) {

	}

	void open() {
		startSimpleGraphRecord(this->out(), this->name());
	}

	void AddVertex(VertexId v) {
		Vertex<VertexId> vertex(v, this->labeler().label(v), this->colorer().GetColour(v));
		recordVertex<VertexId>(this->out(), vertex);
	}

	void AddVertices(const set<VertexId> &elements) {
		map<VertexId, string> colours = this->colorer().GetColours(elements);
		for (auto it = colours.begin(); it != colours.end(); ++ it) {
			Vertex<VertexId> vertex(it->first, this->labeler().label(it->first), it->second);
			recordVertex<VertexId>(this->out(), vertex);
		}
	}

	void AddEdge(EdgeId e, const int length = 0) {
		Edge<VertexId> edge(this->g().EdgeStart(e), this->g().EdgeEnd(e)
				, this->labeler().label(e), this->colorer().GetColour(e), length);
		recordSimpleEdge<VertexId> (this->out(), edge);
	}

	void close() {
		endGraphRecord(this->out());
	}

};

//todo add support of different colors for conjugate vertices
template<class Graph>
class DotPairedGraphPrinter: public GraphPrinter<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef GraphPrinter<Graph> base;
	PairedGraphPrinter<VertexId> paired_printer_;
public:
	DotPairedGraphPrinter(const Graph& g, const GraphLabeler<Graph>& labeler,
			const GraphColorer<Graph>& colorer, const string &name,
			ostream &out = cout) :
		base(g, labeler, colorer, name, out), paired_printer_(name, out) {
	}

	void open() {
		paired_printer_.open();
	}

	void close() {
		paired_printer_.close();
	}

	void AddVertex(VertexId v) {
		VertexId conjugate = this->g().conjugate(v);
		paired_printer_.AddVertex(v, this->labeler().label(v), conjugate
				, this->labeler().label(conjugate), this->colorer().GetColour(v));
	}

	void AddVertices(const set<VertexId> &elements) {
//			map<VertexId, string> colours = this->colorer().GetColours(elements);
//			for (auto it = colours.begin(); it != colours.end(); ++ it) {
//				Vertex<VertexId> vertex(it->first, this->labeler().label(it->first), it->second);
//				recordVertex<VertexId>(this->out(), vertex);
//		}
		set<VertexId> conjugate_vertices;
		for (auto it = elements.begin(); it != elements.end(); ++ it) {
			if (conjugate_vertices.find(*it) == conjugate_vertices.end()){
				AddVertex(*it);
				conjugate_vertices.insert(this->g().conjugate(*it));
			}
		}
	}

	void AddEdge(EdgeId e, const int length = 0) {
		VertexId v1 = this->g().EdgeStart(e);
		VertexId v2 = this->g().EdgeEnd(e);
		paired_printer_.AddEdge(make_pair(v1, this->g().conjugate(v1))
				, make_pair(v2, this->g().conjugate(v2)), this->labeler().label(e)
				, this->colorer().GetColour(e), length);
	}
};

//template<class Graph>
//class PrinterFactory {
//public:
//	virtual auto_ptr<GraphPrinter<Graph>> GetPrinterInstance(
//			const Graph& g, const GraphLabeler<Graph>& labeler,
//			const GraphColorer<Graph>& colorer,
//			const string &graph_name, ostream &os) const = 0;
//	virtual ~PrinterFactory() {
//	}
//};
//
//template<class Graph>
//class DotPrinterFactory {
//public:
//	auto_ptr<GraphPrinter<Graph>> GetPrinterInstance(
//			const Graph& g, const GraphLabeler<Graph>& labeler,
//			const GraphColorer<Graph>& colorer,
//			const string &graph_name, ostream &os) const {
//		return auto_ptr<DotGraphPrinter<Graph>>(new DotGraphPrinter<Graph>(g, labeler, colorer, graph_name, os));
//	}
//};
//
//template<class Graph>
//class PairedDotPrinterFactory {
//public:
//	auto_ptr<GraphPrinter<Graph>> GetPrinterInstance(
//			const Graph& g, const GraphLabeler<Graph>& labeler,
//			const GraphColorer<Graph>& colorer,
//			const string &graph_name, ostream &os) const {
//		return auto_ptr<DotPairedGraphPrinter<Graph>>(new DotPairedGraphPrinter<Graph>(g, labeler, colorer, graph_name, os));
//	}
//};


}
