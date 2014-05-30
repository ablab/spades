//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "standard_base.hpp"
#include "graph_print_utils.hpp"
#include "graph_labeler.hpp"
#include "graph_colorer.hpp"
#include "vertex_linker.hpp"

namespace omnigraph {
namespace visualization {

template<class Graph>
class GraphPrinter {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
//	ostream& os_;
	const Graph &graph_;
protected:
	const GraphLabeler<Graph> &labeler_;
	const GraphColorer<Graph> &colorer_;
	const VertexLinker<Graph> &linker_;

protected:
//	ostream& os() {
//		return os_;
//	}


	const Graph &graph() {
		return graph_;
	}

	template<class GvisVertexId>
	gvis::BaseVertex<GvisVertexId> CreateBaseVertex(GvisVertexId id, VertexId v) {
		return gvis::BaseVertex<GvisVertexId>(id, labeler_.label(v), linker_.GetValue(v), colorer_.GetValue(v));
	}

	template<class GvisVertexId>
	gvis::BaseEdge<GvisVertexId> CreateBaseEdge(GvisVertexId from, GvisVertexId to, EdgeId e){
		return gvis::BaseEdge<GvisVertexId>(from, to, this->labeler_.label(e), this->colorer_.GetValue(e));
	}

	virtual void ManageDrawn(VertexId v, set<VertexId> &visited) {
		visited.insert(v);
	}

public:
	GraphPrinter(const Graph &graph, /*ostream &os,*/
			const GraphLabeler<Graph> &labeler,
			const GraphColorer<Graph> &colorer,
			const VertexLinker<Graph> &linker) :
			/*os_(os), */graph_(graph), labeler_(labeler), colorer_(colorer), linker_(
					linker) {
	}

	virtual void open() = 0;

	virtual void close() = 0;

	virtual void AddVertex(VertexId v1) = 0;

	template<class iter>
	void AddVertices(iter vbegin, iter vend) {
		set<VertexId> drawn;
		for(;vbegin != vend; ++vbegin) {
			if(drawn.count(*vbegin) == 0) {
				AddVertex(*vbegin);
				ManageDrawn(*vbegin, drawn);
			}
		}
	}

	virtual void AddEdge(EdgeId e) = 0;

	template<class iter>
	void AddEdges(iter ebegin, iter eend) {
		for(;ebegin != eend; ++ebegin) {
			AddEdge(*ebegin);
		}
	}

	virtual ~GraphPrinter() {
	}
};

template<typename Graph>
class SingleGraphPrinter : public GraphPrinter<Graph> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	gvis::DotSingleGraphRecorder<size_t> recorder_;

public:
	SingleGraphPrinter(const Graph &graph, ostream &os,
			const GraphLabeler<Graph> &labeler,
			const GraphColorer<Graph> &colorer,
			const VertexLinker<Graph> &linker) : GraphPrinter<Graph>(/*os_, */graph, labeler, colorer, linker), recorder_(os){
	}

	void open() {
		recorder_.startGraphRecord("graph_picture");
	}

	void close() {
		recorder_.endGraphRecord();
	}

	void AddVertex(VertexId v) {
		recorder_.recordVertex(this->CreateBaseVertex((size_t)this->graph().int_id(v), v));
	}

	void AddEdge(EdgeId edge) {
		recorder_.recordEdge(this->CreateBaseEdge((size_t)this->graph().int_id(this->graph().EdgeStart(edge)), (size_t)this->graph().int_id(this->graph().EdgeEnd(edge)), edge));
	}
};

template<typename Graph>
class PairedGraphPrinter : public GraphPrinter<Graph> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	gvis::DotPairedGraphRecorder<size_t> recorder_;

	pair<gvis::BaseVertex<size_t>, gvis::BaseVertex<size_t>> CreateDoubleVertex(VertexId v) {
		gvis::BaseVertex<size_t> u1 = this->CreateBaseVertex((size_t)this->graph().int_id(v), v);
		gvis::BaseVertex<size_t> u2 = this->CreateBaseVertex((size_t)this->graph().int_id(this->graph().conjugate(v)), this->graph().conjugate(v));
		return make_pair(u1, u2);
	}

	pair<size_t, size_t> CreateDoubleVertexId(VertexId v) {
		return make_pair(this->graph().int_id(v), this->graph().int_id(this->graph().conjugate(v)));
	}
protected:
	/*virtual */void ManageDrawn(VertexId v, set<VertexId> &visited) {
		visited.insert(v);
		visited.insert(this->graph().conjugate(v));
	}

public:
	PairedGraphPrinter(const Graph &graph, ostream &os,
			const GraphLabeler<Graph> &labeler,
			const GraphColorer<Graph> &colorer,
			const VertexLinker<Graph> &linker) : GraphPrinter<Graph>(/*os_, */graph, labeler, colorer, linker), recorder_(os) {
	}

	void open() {
		recorder_.startGraphRecord("graph_picture");
	}

	void close() {
		recorder_.endGraphRecord();
	}

	void AddVertex(VertexId v) {
		recorder_.recordVertex(CreateDoubleVertex(v));
	}

	void AddEdge(EdgeId edge) {
		auto vid1 = CreateDoubleVertexId(this->graph().EdgeStart(edge));
		auto vid2 = CreateDoubleVertexId(this->graph().EdgeEnd(edge));
		recorder_.recordEdge(gvis::BaseEdge<pair<size_t, size_t>>(vid1, vid2, this->labeler_.label(edge), this->colorer_.GetValue(edge)));
	}
};

}
}
